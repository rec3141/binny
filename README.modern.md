# binny, modernised front end (branch `modern`)

This fork keeps binny's binning algorithm (Hickl et al., *Briefings in Bioinformatics* 2022) unchanged and replaces
the installation and pre-processing layer:

| upstream (Snakemake workflow) | here |
|---|---|
| Snakemake 7, per-rule conda envs created at run time | one Python package: `binny prepare-db`, `binny run` |
| Prokka (`--metagenome`, `--noanno`) gene calls | Pyrodigal, metagenome mode |
| Mantis, CheckM marker HMMs, E ≤ 1e-3, heuristic domain choice | pyhmmer hmmsearch, same HMMs and E-value, best non-overlapping hits per protein, DUF families dropped |
| bedtools + calcAvgCoverage.pl from BAMs | a depth table (binny format or MetaBAT `jgi_summarize_bam_contig_depths`) |
| barrnap/minced features from Prokka for masking | optional `--mask-gff` (rRNA and CRISPR features) |
| scikit-bio (CLR, multiplicative replacement) | 15 lines of NumPy, same definitions |
| Python 3.8, pandas 1.4, openTSNE 0.6, hdbscan 0.8.30 | Python ≥ 3.10, pandas ≥ 2.1 (tested 3.1), openTSNE ≥ 1.0, hdbscan ≥ 0.8.33, NumPy 2 |

```
conda env create -f envs/binny.yml && pip install .
binny prepare-db checkm_data_2015_01_16/ binny_db/     # unpacked CheckM v1 data
binny run --assembly asm.fa --depth depth.tsv --db binny_db/ --outdir out/ --threads 16
```

Defaults equal `config/config.default.yaml`. The Snakemake workflow (`binny`, `Snakefile`, `workflow/`) is left in
place for comparison but no longer works on this branch (its scripts moved to `src/binny/`).

**Same algorithm as the latest release.** Bioconda's binny 2.2.18 is built from `ohickl/binny` (tag 2.2.18, Aug 2023), a
module-split version of this code. Of the 45 algorithm functions both share, 44 are identical (compared as Python
syntax trees); the 45th differs by a renamed helper.

**Validation on `test/contigs_4bins.fa`** (412 contigs, 4 genomes, 2 depth samples), against binny 2.2.18 from
bioconda run on the same assembly and depth table:

| | binny 2.2.18 (Prokka + Mantis) | this branch (pyrodigal + pyhmmer) |
|---|---|---|
| bins | 4 | 4 |
| contig assignment | — | identical (adjusted Rand index 1.000, all 412 contigs) |
| CheckM marker hits | 8,428 | 8,031 (82% shared per contig) |
| completeness / purity labels | C96-97, P96-98 | C96-100, P97-98 |
| setup | `binny -s`: CheckM download, Mantis + NLTK data, Prokka container; failed without manually installed NLTK data | `binny prepare-db` from CheckM data |
| wall time, 8 threads | 7 min (after setup) | 4.8 min |

The marker search (11,391 CheckM HMMs against every protein) dominates run time in both. A CAMI-scale comparison is
the next check.
