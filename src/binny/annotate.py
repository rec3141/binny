"""Inputs for the binny algorithm without Prokka, Mantis or Snakemake.

- Assembly: contigs of at least `min_length` bp (the Snakemake workflow used `seqkit seq -m 499`).
- Genes: Pyrodigal in metagenome mode (Prokka's `--metagenome` also runs Prodigal in meta mode).
- CheckM markers: pyhmmer hmmsearch of the proteins against the CheckM marker HMMs kept by `binny prepare-db`,
  E-value <= 1e-3 (the Mantis call used `-et 1e-3`). Like Mantis' heuristic domain algorithm, each protein keeps
  its best-scoring hits that do not overlap an already accepted hit by more than `max_overlap` of the shorter
  envelope.
- Annotation GFF in the layout `binny`'s functions read: one CDS line per gene, `ID=<gene>`, plus
  `checkm_marker=<acc,acc>` on marker genes. Optional rRNA / CRISPR features (from any GFF whose rRNA lines carry
  `product=...ribosomal RNA` or whose CRISPR lines carry `rpt_family=`) are appended so binny can mask them from
  k-mer counting, as it did with Prokka's barrnap/minced output.
- Depth: a headerless table `contig<TAB>depth_sample1<TAB>...`. MetaBAT/jgi_summarize_bam_contig_depths tables
  (header `contigName contigLen totalAvgDepth <bam> <bam>-var ...`) are converted.
"""
import gzip
import os
from concurrent.futures import ThreadPoolExecutor


def read_fasta(path):
    op = gzip.open if str(path).endswith('.gz') else open
    name, seq = None, []
    with op(path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, ''.join(seq)


def format_assembly(assembly, out_fa, min_length=499):
    n = kept = 0
    with open(out_fa, 'w') as fo:
        for name, seq in read_fasta(assembly):
            n += 1
            if len(seq) >= min_length:
                kept += 1
                fo.write(f'>{name}\n')
                for i in range(0, len(seq), 80):
                    fo.write(seq[i:i + 80] + '\n')
    return n, kept


def convert_depth(depth_file, out_file, contigs=None):
    """write binny's headerless depth table; returns the number of samples"""
    op = gzip.open if str(depth_file).endswith('.gz') else open
    with op(depth_file, 'rt') as f:
        lines = [l.rstrip('\n').split('\t') for l in f if l.strip()]
    header = lines[0]
    if header[0] == 'contigName' and len(header) > 2 and header[2] == 'totalAvgDepth':
        cols = [i for i in range(3, len(header)) if not header[i].endswith('-var')]
        rows = [[r[0]] + [r[i] for i in cols] for r in lines[1:]]
    else:
        try:
            [float(x) for x in header[1:]]
            rows = lines
        except ValueError:          # some other header line: drop it
            rows = lines[1:]
    keep = set(contigs) if contigs is not None else None
    n = 0
    with open(out_file, 'w') as fo:
        for r in rows:
            if keep is None or r[0] in keep:
                fo.write('\t'.join(r) + '\n'); n += 1
    return len(rows[0]) - 1, n


def call_genes(contigs, threads=1):
    """Pyrodigal, metagenome mode. Returns list of (gene_id, contig, start, end, strand, protein)."""
    import pyrodigal
    finder = pyrodigal.GeneFinder(meta=True)

    def one(item):
        name, seq = item
        out = []
        for k, g in enumerate(finder.find_genes(seq.encode())):
            out.append((f'{name}_{k + 1}', name, g.begin, g.end, '+' if g.strand == 1 else '-',
                        g.translate().rstrip('*')))
        return out
    genes = []
    with ThreadPoolExecutor(max_workers=max(1, threads)) as ex:     # pyrodigal releases the GIL
        for res in ex.map(one, contigs, chunksize=64):
            genes.extend(res)
    return genes


def find_markers(genes, hmm_files, threads=1, evalue=1e-3, max_overlap=0.1):
    """pyhmmer hmmsearch; returns {gene_id: [marker accession, ...]} with versions stripped"""
    import pyhmmer
    alphabet = pyhmmer.easel.Alphabet.amino()
    seqs = pyhmmer.easel.TextSequenceBlock(
        pyhmmer.easel.TextSequence(name=g[0].encode(), sequence=g[5]) for g in genes if g[5]).digitize(alphabet)
    hits = {}
    for hf in hmm_files:
        with pyhmmer.plan7.HMMFile(hf) as fh:
            hmms = list(fh)
        for top in pyhmmer.hmmsearch(hmms, seqs, cpus=threads, E=evalue):
            q = top.query
            name_q = q.name.decode() if isinstance(q.name, bytes) else q.name
            if name_q.startswith('DUF'):       # binny ignored DUF annotations (parse_mantis_cons_annot)
                continue
            acc = (q.accession or q.name)
            acc = (acc.decode() if isinstance(acc, bytes) else acc).split('.')[0]
            for h in top:
                if not h.included:
                    continue
                d = h.best_domain
                name = h.name.decode() if isinstance(h.name, bytes) else h.name
                hits.setdefault(name, []).append((h.score, d.env_from, d.env_to, acc))
    out = {}
    for gene, hs in hits.items():
        hs.sort(key=lambda x: -x[0])
        acc_list, spans = [], []
        for score, a, b, acc in hs:
            ok = True
            for c, e in spans:
                ov = min(b, e) - max(a, c) + 1
                if ov > max_overlap * min(b - a + 1, e - c + 1):
                    ok = False; break
            if ok and acc not in acc_list:
                acc_list.append(acc); spans.append((a, b))
        out[gene] = acc_list
    return out


def write_gff(genes, markers, out_gff, mask_gff=None):
    with open(out_gff, 'w') as fo:
        for gid, contig, start, end, strand, _ in genes:
            attr = f'ID={gid}'
            if markers.get(gid):
                attr += ';checkm_marker=' + ','.join(markers[gid])
            fo.write(f'{contig}\tpyrodigal\tCDS\t{start}\t{end}\t.\t{strand}\t0\t{attr}\n')
        if mask_gff:
            with open(mask_gff) as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    p = line.rstrip('\n').split('\t')
                    if len(p) == 9 and (('product=' in p[8] and 'ribosomal RNA' in p[8]) or 'rpt_family=' in p[8]):
                        fo.write(line if line.endswith('\n') else line + '\n')
