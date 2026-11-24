"""binny command line (modernised front end; the binning algorithm is unchanged).

  binny prepare-db CHECKM_DATA_DIR DB_DIR
      CHECKM_DATA_DIR: unpacked CheckM v1 data (checkm_data_2015_01_16: taxon_marker_sets.tsv, pfam/tigrfam2pfam.tsv,
      hmms/checkm.hmm).
  binny run --assembly ASM.fa --depth DEPTH.tsv --db DB_DIR --outdir OUT [--threads N] [--mask-gff rRNA_CRISPR.gff]
      DEPTH.tsv: contig<TAB>depth per sample (no header), or a MetaBAT jgi_summarize_bam_contig_depths table.
"""
import argparse
import os
import sys
import time


def _bool(x):
    return str(x).lower() in ('1', 'true', 'yes', 'on')


def main(argv=None):
    ap = argparse.ArgumentParser(prog='binny', description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('prepare-db', help='build the CheckM marker database')
    p.add_argument('checkm_data'); p.add_argument('db')
    r = sub.add_parser('run', help='bin an assembly')
    r.add_argument('--assembly', required=True)
    r.add_argument('--depth', required=True)
    r.add_argument('--db', required=True)
    r.add_argument('--outdir', required=True)
    r.add_argument('--sample', default='sample')
    r.add_argument('--threads', type=int, default=1)
    r.add_argument('--mask-gff', help='GFF with rRNA (product=...ribosomal RNA) and CRISPR (rpt_family=) features')
    r.add_argument('--kmers', default='2,3,4')
    r.add_argument('--mask-disruptive-sequences', type=_bool, default=True)
    r.add_argument('--extract-scmags', type=_bool, default=True)
    r.add_argument('--coassembly-mode', choices=['auto', 'on', 'off'], default='auto')
    r.add_argument('--nx-value', type=int, default=90)
    r.add_argument('--min-cont-length-cutoff', type=int, default=2250)
    r.add_argument('--max-cont-length-cutoff', type=int, default=2250)
    r.add_argument('--min-cont-length-cutoff-marker', type=int, default=2250)
    r.add_argument('--max-cont-length-cutoff-marker', type=int, default=2250)
    r.add_argument('--max-n-contigs', type=float, default=5.0e5)
    r.add_argument('--max-marker-lineage-depth-lvl', type=int, default=2)
    r.add_argument('--distance-metric', default='manhattan')
    r.add_argument('--max-iterations', type=int, default=50)
    r.add_argument('--hdbscan-epsilon-range', default='0.250,0.000')
    r.add_argument('--hdbscan-min-samples-range', default='1,5,10')
    r.add_argument('--include-depth-initial', type=_bool, default=False)
    r.add_argument('--include-depth-main', type=_bool, default=False)
    r.add_argument('--min-completeness', type=float, default=72.5)
    r.add_argument('--start-completeness', type=float, default=92.5)
    r.add_argument('--purity', type=float, default=95)
    r.add_argument('--write-contig-data', type=_bool, default=True)
    a = ap.parse_args(argv)

    if a.cmd == 'prepare-db':
        from .checkm_db import prepare_db
        prepare_db(a.checkm_data, a.db)
        print(f'binny database written to {a.db}')
        return 0

    from . import annotate
    out = os.path.abspath(a.outdir)
    inter = os.path.join(out, 'intermediary')
    os.makedirs(inter, exist_ok=True); os.makedirs(os.path.join(out, 'logs'), exist_ok=True)
    asm = os.path.join(inter, 'assembly.formatted.fa')
    depth = os.path.join(inter, 'assembly.contig_depth.txt')
    gff = os.path.join(inter, 'annotation_CDS_RNA_hmms_checkm.gff')
    t0 = time.time()

    def say(msg):
        print(f'[{time.time() - t0:7.0f}s] {msg}', file=sys.stderr, flush=True)
    # the Snakemake workflow kept contigs >= min_cont_length_cutoff_marker - 1 (seqkit -m 499 with the 500 default)
    n, kept = annotate.format_assembly(a.assembly, asm, min_length=499)
    say(f'assembly: {kept} of {n} contigs >= 499 bp')
    contigs = [c for c, _ in annotate.read_fasta(asm)]
    ns, nd = annotate.convert_depth(a.depth, depth, contigs)
    say(f'depth: {ns} sample(s), {nd} contigs')
    if not os.path.exists(gff):
        genes = annotate.call_genes(list(annotate.read_fasta(asm)), a.threads)
        say(f'genes: {len(genes)}')
        hmms = [os.path.join(a.db, 'hmms', 'checkm_pf', 'checkm_filtered_pf.hmm'),
                os.path.join(a.db, 'hmms', 'checkm_tf', 'checkm_filtered_tf.hmm')]
        markers = annotate.find_markers(genes, hmms, a.threads)
        say(f'CheckM markers on {len(markers)} genes')
        annotate.write_gff(genes, markers, gff + '.tmp', a.mask_gff)
        os.replace(gff + '.tmp', gff)
    from .pipeline import run
    os.chdir(out)          # the algorithm writes some intermediary files relative to the working directory
    bins = run(out, a.sample, depth, asm, gff,
               os.path.join(a.db, 'pfam', 'tigrfam2pfam.tsv'),
               os.path.join(a.db, 'taxon_marker_sets_lineage_sorted.tsv'),
               threads=a.threads, min_purity=a.purity, min_completeness=a.min_completeness,
               starting_completeness=a.start_completeness, kmers=a.kmers,
               mask_disruptive_sequences=a.mask_disruptive_sequences, extract_scmags=a.extract_scmags,
               coassembly_mode=a.coassembly_mode, nx_val=a.nx_value,
               min_contig_length=a.min_cont_length_cutoff, max_contig_length=a.max_cont_length_cutoff,
               min_contig_length_marker=a.min_cont_length_cutoff_marker,
               max_contig_length_marker=a.max_cont_length_cutoff_marker, max_contig_threshold=a.max_n_contigs,
               max_marker_lineage_depth_lvl=a.max_marker_lineage_depth_lvl, max_embedding_tries=a.max_iterations,
               include_depth_initial=a.include_depth_initial, include_depth_main=a.include_depth_main,
               hdbscan_epsilon_range=a.hdbscan_epsilon_range.split(','),
               hdbscan_min_samples_range=a.hdbscan_min_samples_range.split(','),
               dist_metric=a.distance_metric, write_contig_data=a.write_contig_data,
               log=os.path.join(out, 'logs', 'binny.log'))
    say(f'done: {len(bins)} bins in {os.path.join(out, "bins")}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
