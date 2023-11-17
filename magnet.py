import os
import argparse
import pathlib
import sys
import subprocess
import math
import warnings
from io import StringIO
from collections import defaultdict
from multiprocessing import Pool, Manager
from itertools import repeat
warnings.filterwarnings("ignore")

import pandas as pd
from ete3 import NCBITaxa
from Bio import SeqIO

from utils.reference_finder import prepare_reference_genomes
from utils.alignment import run_minimap2, run_bwa, sort_samfile, samtools_calculate_coverage
from utils.summary import alignment_summary, merge_reference_fasta, call_present_absent
from utils.ani import samtools_merged_consensus, ani_summary
from utils.input_parsing import parsing_input_f, filter_input_df, get_seq2assembly_dict

def main():
    parser = argparse.ArgumentParser(description="Universal Taxonomic Classification Verifier.")

    parser.add_argument("-c", "--classification", type=pathlib.Path, required=True, help="Path to the Taxonomic Classification Report. Accepting csv/tsv file format, other text formats are treated as tsv.")
    parser.add_argument("-i", "--fastq", type=pathlib.Path, required=True, help="Path to the first fastq file.")
    parser.add_argument("-I", "--fastq2", type=pathlib.Path, required=False, help="Path to the second fastq file for paired-end reads.")
    parser.add_argument("-m", "--mode", type=str, required=False, choices=['ont', 'illumina'], help="Modes for different sequencing platforms [ont, illumina]. Default:[ont]",  default='ont')
    parser.add_argument("-o", "--output", type=pathlib.Path, required=True, help="Path to the output directory.")
    parser.add_argument("-t", "--taxid-idx", type=int, required=False, help="The column index (0-based) of the taxids. Default:[0]", default=0)
    parser.add_argument("-a", "--abundance-idx", type=int, required=False, help="The column index (0-based) of the abundance. Default:[None]")
    parser.add_argument("--min-abundance", type=float, required=False, help="Minimum abundance (0-1) for pre-filtering, exclude taxa below the threshold.", default=0)
    parser.add_argument("--min-mapq", type=int, required=False, help="Minimum MAPQ for primary alignments. Default:[20]", default=20)
    parser.add_argument("--min-covscore", type=float, required=False, help="Minimum Coverage Score for supplementary alignments. Default:[0.7]", default=0.7)
    parser.add_argument("--threads", type=int, required=False, help="Number of threads for Multi-threading. Default:[1]", default=1)
    parser.add_argument("--include-mag", action='store_true', required=False, help="Include metagenomic assemble genomes. Default:[off]")
    parser.set_defaults(include_mag=False)
    parser.add_argument("--subspecies", action='store_true', required=False, help="Verify taxonomic classification at subspecies rank. Default:[off]")
    parser.set_defaults(subspecies=False)

    args = parser.parse_args()
    
    input_tsv = args.classification
    input_fastq = args.fastq
    input_fastq2 = args.fastq2
    mode = args.mode
    working_directory = args.output

    taxid_col_idx = args.taxid_idx
    abundance_col_idx = args.abundance_idx
    min_abundance = args.min_abundance
    min_mapq = args.min_mapq
    min_coverage_score = args.min_covscore
    threads = args.threads

    if args.include_mag:
        mag_flag = 'all'
    else:
        mag_flag = 'exclude'

    call_subspecies = args.subspecies

    sep = '\t'
    if str(input_tsv)[-3:] == 'csv':
        sep = ','
        
    ncbi_taxa_db = NCBITaxa()
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)
        
    input_df, min_abundance = parsing_input_f(input_tsv, sep, taxid_col_idx, abundance_col_idx, min_abundance)
    # make valid_kingdom a variable?
    valid_taxids = filter_input_df(input_df, min_abundance, ncbi_taxa_db, valid_kingdom={2, 4751, 2157, 10239}, ret_subspecies=call_subspecies)
    
    # can be parallelized
    reference_metadata = prepare_reference_genomes(valid_taxids, working_directory, ncbi_taxa_db, mag_flag=mag_flag)
    downloaded_assemblies = reference_metadata[reference_metadata['Downloaded']]
    seq2assembly_dict = get_seq2assembly_dict(working_directory, downloaded_assemblies)
    reference_fasta = merge_reference_fasta(list(downloaded_assemblies['Assembly Accession ID']), working_directory)
    
    if mode == 'ont':
        aligner_output = run_minimap2(input_fastq, reference_fasta, 'merged', working_directory, threads=threads)
    if mode == 'illumina':
        aligner_output = run_bwa(input_fastq, input_fastq2, reference_fasta, 'merged', working_directory, threads=threads)
        
    sort_samfile('merged', aligner_output, working_directory, min_mapq=0, threads=threads)
    
    coverage_files = os.path.join(working_directory, "coverage_files")
    if not os.path.exists(coverage_files):
        os.mkdir(coverage_files)
        
    pool = Pool(processes=threads)
    pool.starmap(samtools_calculate_coverage, zip(repeat(working_directory), [True, False]))
    pool.close()
    pool.join()
    
    downloaded_assemblies = alignment_summary(downloaded_assemblies,
                                          working_directory,
                                          seq2assembly_dict,
                                          include_supp=True)
    
    downloaded_assemblies = alignment_summary(downloaded_assemblies,
                                          working_directory,
                                          seq2assembly_dict,
                                          include_supp=False)
    
    consensus_record_dict = samtools_merged_consensus(working_directory, threads)
    downloaded_assemblies = ani_summary(downloaded_assemblies, consensus_record_dict, working_directory, threads)
    downloaded_assemblies = call_present_absent(downloaded_assemblies, min_coverage_score)
    downloaded_assemblies.sort_values(['Primary Score'], ascending=False).to_csv(os.path.join(working_directory, 
                                                                                          'magnet_results.csv'), 
                                                                             index=False)
    
if __name__ == "__main__":
    main()
