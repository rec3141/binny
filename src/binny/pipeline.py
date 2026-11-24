"""binny main loop as a function (formerly workflow/scripts/binny_main.py, a Snakemake script).

Algorithm unchanged; parameters come from the CLI (binny.cli) instead of the Snakemake object, and the annotation
GFF already carries the CheckM markers (binny.annotate), so the Mantis-to-GFF merge step is gone.
Original author: Oskar Hickl (2021).
"""
import glob
import logging
import os

import numpy as np
import pandas as pd

from .functions import *  # noqa: F401,F403  (the original script did the same via sys.path)


def run(binny_out, sample, mg_depth_file, assembly, annot_file, tigrfam2pfam_file, taxon_marker_set_file,
        threads=1, min_purity=95.0, min_completeness=72.5, starting_completeness=92.5, kmers='2,3,4',
        mask_disruptive_sequences=True, extract_scmags=True, coassembly_mode='auto', nx_val=90,
        min_contig_length=2250, max_contig_length=2250, min_contig_length_marker=2250,
        max_contig_length_marker=2250, max_contig_threshold=5.0e5, max_marker_lineage_depth_lvl=2,
        max_embedding_tries=50, include_depth_initial=False, include_depth_main=False,
        hdbscan_epsilon_range=(0.250, 0.000), hdbscan_min_samples_range=(1, 5, 10), dist_metric='manhattan',
        write_contig_data=True, log=None):
    hdbscan_epsilon_range = [float(e) for e in hdbscan_epsilon_range]
    hdbscan_min_samples_range = [int(m) for m in hdbscan_min_samples_range]
    n_dim = 2
    # To achieve reproducible results with HDBSCAN and ensure same seed, because other tools that accept seed
    # arguments, might mess with the global numpy seed
    np.random.seed(0)
    logging.basicConfig(filename=log, level=logging.INFO, format='%(asctime)s - %(message)s',
                        datefmt='%d/%m/%Y %I:%M:%S %p', filemode='w', force=True)
    logging.getLogger('numba').setLevel(logging.WARNING)
    logging.info('Starting Binny run for sample {0}.'.format(sample))

    # Check if bin dir empty
    if glob.glob(os.path.join(binny_out, "bins/*.fasta")):
        logging.error('Bin dir contains fasta files. Move or delete. Exiting.')
        raise Exception

    # Load TIGRFAMs to PFAMs conversion table.
    tigrfam2pfam_data = tigrfam2pfam_dict(tigrfam2pfam_file)

    # Merge Prokka gff with marker set data, load annotation df, and load assembly.
    annot_df, annot_dict = gff2ess_gene_df(annot_file, target_attribute='checkm_marker', get_dict=True)
    assembly_dict = load_fasta(assembly)

    # Build marker set graph db.
    taxon_marker_sets = load_checkm_markers(taxon_marker_set_file)

    # Look for complete genomes on single contigs
    all_good_bins = {}
    if extract_scmags:
        logging.info('Looking for single contig bins.')
        single_contig_bins = get_single_contig_bins(annot_df, all_good_bins, n_dim, taxon_marker_sets, tigrfam2pfam_data,
                                                    threads, max_marker_lineage_depth_lvl=max_marker_lineage_depth_lvl)
        logging.info('Found {0} single contig bins.'.format(len(single_contig_bins)))
    else:
        single_contig_bins = []

    # logging.info('Getting assembly dict without scMAG.'.format(len(single_contig_bins)))
    # assembly_dict_wo_scmags = {key: val for key, val in assembly_dict.items() if key in
    #                            set(assembly_dict.keys()).difference(set(single_contig_bins))}

    logging.info(f'Calculating N{nx_val}'.format(len(single_contig_bins)))
    nx = calc_assembly_nx(assembly_dict, single_contig_bins, nx_val)
    nx2 = calc_assembly_nx(assembly_dict, [], nx_val)
    logging.info(f'N{nx_val} is {nx}, with scMAGs would be {nx2}.'.format(len(single_contig_bins)))
    min_contig_length = min(max(nx, min_contig_length), max_contig_length)
    min_contig_length_marker = min(max(int(nx / 3), min_contig_length_marker), max_contig_length_marker)

    # Load assembly and mask rRNAs and CRISPR arrays
    # contig_list = [[contig] + [seq] for contig, seq in assembly_dict.items() if (len(seq) >= min_contig_length
    #                                                                              or (annot_dict.get(contig)
    #                                                                              and len(seq) >= min_contig_length_marker))
    #                                                                             and contig not in single_contig_bins]
    contig_list = [[contig] + [seq] for contig, seq in assembly_dict.items() if contig not in single_contig_bins
                                                                                and len(seq) >= 500]

    # logging.info('{0} contigs match length threshold of {1}bp or contain marker genes and'
    #              ' have a size of at least {2}bp'.format(len(contig_list), min_contig_length, min_contig_length_marker))

    contig_rrna_crispr_region_dict = gff2low_comp_feature_dict(annot_file)
    if mask_disruptive_sequences:
        logging.info('Masking potentially disruptive sequences from k-mer counting.')
        mask_rep_featrues(contig_rrna_crispr_region_dict, contig_list)
    else:
        logging.info('Not masking potentially disruptive sequences from k-mer counting.')

    # Get length normalized k-mer frequencies.
    kmer_sizes = [int(kmer) for kmer in kmers.split(',')]

    start = timer()
    kfreq_array = get_contig_kmer_matrix(contig_list, kmer_sizes, threads)
    end = timer()
    logging.info('K-mer frequency matrix created in {0}s.'.format(int(end - start)))

    # Make array, removing fully masked sequences with no counts and standardize k-mer freq data
    x = np.array([c_kfreq[1:] for c_kfreq in kfreq_array[1:] if not sum(c_kfreq[1:]) == 0])
    x_contigs = [c_kfreq[0] for c_kfreq in kfreq_array[1:] if not sum(c_kfreq[1:]) == 0]

    main_contig_data_dict = {cont: seq for cont, seq in zip(x_contigs, x)}

    # Load depth data
    depth_dict = load_depth_dict(mg_depth_file)
    n_depth_samples = list(depth_dict.values())[0].shape[0]
    logging.debug(f'list(depth_dict.values())[0].shape[0]: {n_depth_samples}')

    if coassembly_mode == 'on' or (coassembly_mode == 'auto' and n_depth_samples > 1):
        min_contig_length_marker = 500
        if n_depth_samples >= 3:
            include_depth_main = True
        logging.info('Using coassembly mode.')
    else:
        logging.info('Using single sample mode.')

    # Run iterative dimension reduction, manifold learning, cluster detection and assessment.
    all_good_bins, contig_data_df_org, min_purity = iterative_embedding(x_contigs, depth_dict, all_good_bins, starting_completeness,
                                                            min_purity, min_completeness, threads, n_dim, annot_file,
                                                            mg_depth_file, single_contig_bins, taxon_marker_sets,
                                                            tigrfam2pfam_data, main_contig_data_dict, assembly_dict,
                                                            max_contig_threshold, min_contig_length_marker,
                                                            include_depth_initial, max_embedding_tries,
                                                            include_depth_main, hdbscan_epsilon_range,
                                                            hdbscan_min_samples_range, dist_metric,
                                                            contigs2clusters_out_path=os.path.join(binny_out, 'intermediary'),
                                                            max_marker_lineage_depth_lvl=max_marker_lineage_depth_lvl,
                                                            coassembly_mode=coassembly_mode)

    all_contigs = []
    for bin in all_good_bins:
        all_contigs.extend(all_good_bins[bin]['contigs'])
    if len(all_contigs) != len(set(all_contigs)):
        logging.warning('WARNING: {0} duplicate contigs in bins found!'.format(len(all_contigs) - len(set(all_contigs))))

    # Write bin fastas.
    write_bins(all_good_bins, assembly, min_comp=int(min_completeness), min_pur=int(min_purity),
               bin_dir=os.path.join(binny_out, 'bins'))

    bin_dict = {contig: bin for bin, bin_data in all_good_bins.items() for contig in bin_data['contigs']}

    all_cont_data_dict = {}

    for contig, k_freq in main_contig_data_dict.items():
        all_cont_data_dict[contig] = {'bin': bin_dict.get(contig, 'N'),
                                      'k-mer_freqs': ';'.join([str(k) for k in list(k_freq)]),
                                      'depths': ';'.join([str(d) for d in list(depth_dict.get(contig))])}

    all_cont_data_dict.update({contig: {'bin': contig,
                                      'k-mer_freqs': '',
                                      'depths': ';'.join([str(d) for d in list(depth_dict.get(contig))])}
                               for contig in single_contig_bins})

    contig_data_df = pd.DataFrame.from_dict(all_cont_data_dict, orient='index', columns=['bin', 'k-mer_freqs', 'depths'])
    # Downcast df
    downcast_pd_df(contig_data_df)

    if write_contig_data:
        logging.info('Writing contig data to file.')
        compression_opts = dict(method='gzip')
        contig_data_df.to_csv(os.path.join(binny_out, 'contig_data.tsv.gz'), header=True, index=True, index_label='contig',
                              chunksize=250000, compression=compression_opts, sep='\t')

    open(os.path.join(binny_out, 'binny.done'), 'w').close()
    return all_good_bins

    logging.info('Run finished.')
