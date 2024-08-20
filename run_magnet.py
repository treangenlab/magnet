import os
import argparse
import pathlib
import sys
import subprocess
import math
import warnings
from collections import defaultdict
from multiprocessing import Pool, Manager
from itertools import repeat
warnings.filterwarnings("ignore")

import pandas as pd
from ete3 import NCBITaxa
from Bio import SeqIO

from utils.reference_finder import prepare_reference_genomes
from utils.alignment import run_minimap2, run_bwa, sort_samfile
from utils.summary import merge_reference_fasta, call_present_absent
from utils.ani import ani_summary
from utils.input_parsing import parsing_input_f, filter_input_df, get_seq2assembly_dict

import numpy as np
from sklearn.cluster import AgglomerativeClustering
from ast import literal_eval
    
def samtools_calculate_coverage(output_dir, include_supp=False):
    coverage_files = os.path.join(output_dir, "coverage_files")
    bam_files = os.path.join(output_dir, "bam_files")

    command = ["samtools",
               "coverage", "--no-header",
               os.path.join(bam_files, f"merged.sorted.bam")]
    
    if include_supp:
        # samtools coverage --ff UNMAP,QCFAIL,DUP -q 0 merged.sorted.bam > secondary_coverage.tsv
        coverage_file = os.path.join(coverage_files, f"secondary_coverage.tsv")
        command += ['--ff', 'UNMAP,QCFAIL,DUP', '-q', str(1)]
    else:
        # samtools coverage merged.sorted.bam -q 20 > primary_coverage.tsv
        coverage_file = os.path.join(coverage_files, f"primary_coverage.tsv")
        command += ['-q', str(20)]
        # command += ['--ff', 'UNMAP,QCFAIL,DUP', '-q', str(20)]

    subprocess.run(command,
                    check=True,
                    stdout=open(coverage_file, "w"))
    
def parse_fastani_line(line):
    ret_list = []
    for item in line.strip().split('\t')[1:]:
        if item == 'NA':
            ret_list.append(float(0))
        else:
            ret_list.append(float(item))
    return ret_list

def find_representative_genome(fastani_path, fastani_assemblies, downloaded_assemblies):
    fastani_result = os.path.join(fastani_path, f'pairwise_ani.matrix')
    #fastani_assemblies = downloaded_assemblies[downloaded_assemblies['Genus Taxid'] == genus_taxid]['Assembly Accession ID'].values
    

    ani_matrix = []
    with open(fastani_result, 'r') as fastani_out_f:
        lines = fastani_out_f.readlines()
        num_seqs = int(lines[0].strip())
        
        for idx, line in enumerate(lines[1:]):
            ani_matrix.append(parse_fastani_line(line)+list(np.ones(num_seqs-idx)*100))

    dist_nparray = 1-np.array(ani_matrix)/100
    dist_df = pd.DataFrame(dist_nparray.T + dist_nparray,
                 columns=fastani_assemblies,
                 index=fastani_assemblies)

    model = AgglomerativeClustering(affinity='precomputed', n_clusters=None, compute_full_tree=True,
                                    linkage='complete', 
                                    distance_threshold=0.05).fit(dist_df)
    cluster_df = dist_df.copy()
    #print(cluster_df)
    cluster_df['Cluster Label'] = model.labels_
    
    representative_genomes = defaultdict(list)
    member2representative = dict()
    for cluster_idx in cluster_df['Cluster Label'].unique():
        cluster_members = cluster_df[cluster_df['Cluster Label'] == cluster_idx].index.to_list()
        selected_df = downloaded_assemblies[downloaded_assemblies['Assembly Accession ID'].isin(cluster_members)].copy()
        if selected_df[selected_df['Assembly Level'] == 'Complete Genome'].shape[0] == 1:
            representative_genome = selected_df[selected_df['Assembly Level'] == 'Complete Genome']['Assembly Accession ID'].values[0]
        elif selected_df[selected_df['Assembly Level'] == 'Complete Genome'].shape[0] > 1:
            complete_genomes = list(selected_df[selected_df['Assembly Level'] == 'Complete Genome']['Assembly Accession ID'].values)
            representative_genome = dist_df.loc[complete_genomes].sum().idxmin()
        else:
            representative_genome = dist_df.loc[cluster_members].sum().idxmin()
        
        representative_genomes[representative_genome] = cluster_members
        for member in cluster_members:
            member2representative[member] = representative_genome
            
    return representative_genomes, member2representative

def samtools_merged_consensus(output_directory, threads):
    merged_bam = os.path.join(output_directory, 'bam_files', 'merged.sorted.bam')
    subprocess.run(['samtools', 'consensus', 
                    '--show-ins', 'no', 
                    '--show-del', 'yes',
                    '--min-MQ', str(20),
                    '-a',
                    '--mode', "simple",
                    '--threads', str(threads),
                    merged_bam, 
                    '-o', os.path.join(output_directory, 'merged_consensus.fasta')],
                  check=True)
    
    consensus_record_dict = SeqIO.to_dict(SeqIO.parse(os.path.join(output_directory, 'merged_consensus.fasta'), "fasta"))
    return consensus_record_dict

def alignment_summary(downloaded_assemblies, output_directory, seq2assembly_dict, include_supp=True):
    if include_supp:
        coverage_file_name = 'secondary_coverage.tsv'
        column_prefix = 'Secondary'
    else:
        coverage_file_name = 'primary_coverage.tsv'
        column_prefix = 'Primary'
    
    columns = [f'{column_prefix} Breadth',
               f'{column_prefix} Expected',
               f'{column_prefix} Score',
               f'{column_prefix} Depth']
    
    coverage_df = pd.read_csv(os.path.join(output_directory,
                                           'coverage_files',
                                           coverage_file_name),
                              sep='\t',
                              header=None,
                              names=['rname','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq'])
    
    taxa_records = defaultdict(lambda: defaultdict(int))
    for idx, row in coverage_df.iterrows():
        taxa_reference = seq2assembly_dict[row['rname']]
        taxa_records[taxa_reference]['genome_length'] += row['endpos']
        taxa_records[taxa_reference]['reads_mapped'] += row['numreads']
        taxa_records[taxa_reference]['genome_totol_count'] += int(row['meandepth'] * row['endpos'])
        taxa_records[taxa_reference]['covbases'] += row['covbases']
    
    breadth_coverage_list = []
    depth_coverage_list = []
    expected_breadth_coverage_list = []
    coverage_score = []
    for assembly_id in downloaded_assemblies['Assembly Accession ID']:
        breadth_coverage, depth_coverage, expected_breadth_coverage = calculate_depth(assembly_id, taxa_records)
        breadth_coverage_list.append(breadth_coverage)
        depth_coverage_list.append(depth_coverage)
        expected_breadth_coverage_list.append(expected_breadth_coverage)
        if expected_breadth_coverage != 0:
            coverage_score.append(min(breadth_coverage/expected_breadth_coverage, 1))
        else:
            coverage_score.append(0)
        
    downloaded_assemblies[columns[0]] = breadth_coverage_list
    downloaded_assemblies[columns[1]] = expected_breadth_coverage_list
    downloaded_assemblies[columns[2]] = coverage_score
    downloaded_assemblies[columns[3]] = depth_coverage_list
    
    downloaded_assemblies.to_csv(os.path.join(output_directory, 'alignment.csv'), index=False)
    
    return downloaded_assemblies

def run_magnet(cmd_args):
    args = parser.parse_args(cmd_args)

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
    valid_kingdom_str = args.kingdom

    valid_kingdom = set()
    for i in valid_kingdom_str.split(','):
        valid_kingdom.add(int(i))

    if args.include_mag:
        mag_flag = 'all'
    else:
        mag_flag = 'exclude'

    accession_flag = args.accession
    call_subspecies = args.subspecies

    sep = '\t'
    if str(input_tsv)[-3:] == 'csv':
        sep = ','

    ncbi_taxa_db = NCBITaxa()
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    if accession_flag:
        abundance_col_idx = None
        min_abundance = 0
        input_df, min_abundance = parsing_input_f(input_tsv, sep, taxid_col_idx, abundance_col_idx, min_abundance)
        valid_taxids = list(input_df['tax_id'].values)
    else:
        input_df, min_abundance = parsing_input_f(input_tsv, sep, taxid_col_idx, abundance_col_idx, min_abundance)
        # make valid_kingdom a variable?
        valid_taxids = filter_input_df(input_df, min_abundance, ncbi_taxa_db, valid_kingdom=valid_kingdom, ret_subspecies=call_subspecies)

    reference_metadata = prepare_reference_genomes(valid_taxids, working_directory, ncbi_taxa_db, accession_flag=accession_flag, mag_flag=mag_flag)
    downloaded_assemblies = reference_metadata[reference_metadata['Downloaded']]

    reference_genome_path = os.path.join(working_directory, 'reference_genomes')

    fastani_path = os.path.join(working_directory, 'fastANI')
    if not os.path.exists(fastani_path):
        os.mkdir(fastani_path)

    fastani_assemblies = downloaded_assemblies['Assembly Accession ID'].values
    assemblie_list_f = os.path.join(fastani_path, f"assemblie_list.txt")
    with open(assemblie_list_f, "w") as rl_f:
        for accession in fastani_assemblies:
            reference_genome = os.path.join(reference_genome_path, f'{accession}.fasta')
            rl_f.write(f"{reference_genome}\n")

    subprocess.run(['fastANI',
                    '--rl', assemblie_list_f,
                    '--ql', assemblie_list_f,
                    "--matrix",
                    '--threads', str(threads),
                    '-o', os.path.join(fastani_path, f'pairwise_ani')],
                   check=True,
                   stdout=open(os.path.join(fastani_path, "fastani.log"), "a"),
                   stderr=open(os.path.join(fastani_path, "fastani.err"), "a"))

    representative_genomes, member2representative = find_representative_genome(fastani_path, fastani_assemblies, downloaded_assemblies)

    representative_labels = []
    cluster_members = []
    for idx, row in downloaded_assemblies.iterrows():
        accession = row['Assembly Accession ID']
        if accession in representative_genomes.keys():
            representative_labels.append(True)
            cluster_members.append(','.join(representative_genomes[accession]))
        else:
            representative_labels.append(False)
            cluster_members.append(','.join(representative_genomes[member2representative[accession]]))
    downloaded_assemblies['Cluster Representative'] = representative_labels
    downloaded_assemblies['Cluster Members'] = cluster_members

    representative_df = downloaded_assemblies[downloaded_assemblies['Cluster Representative']]

    seq2assembly_dict = get_seq2assembly_dict(working_directory, representative_df)
    reference_fasta = merge_reference_fasta(list(representative_df['Assembly Accession ID']), working_directory)

    if mode == 'ont':
        aligner_output = run_minimap2(input_fastq, reference_fasta, 'merged', working_directory, threads=threads)
    if mode == 'illumina':
        aligner_output = run_bowtie2(input_fastq, input_fastq2, reference_fasta, 'merged', working_directory, threads=threads)
    sort_samfile('merged', aligner_output, working_directory, min_mapq=0, threads=threads)

    coverage_files = os.path.join(working_directory, "coverage_files")
    if not os.path.exists(coverage_files):
        os.mkdir(coverage_files)

    pool = Pool(processes=threads)
    pool.starmap(samtools_calculate_coverage, zip(repeat(working_directory), [True, False]))
    pool.close()
    pool.join()


    representative_df = alignment_summary(representative_df,
                                          working_directory,
                                          seq2assembly_dict,
                                          include_supp=True)

    representative_df = alignment_summary(representative_df,
                                          working_directory,
                                          seq2assembly_dict,
                                          include_supp=False)

    consensus_record_dict = samtools_merged_consensus(working_directory, threads)
    representative_df = ani_summary(representative_df, consensus_record_dict, working_directory, threads)
    representative_df = call_present_absent(representative_df, min_coverage_score)
    representative_df.sort_values(['Primary Score'], ascending=False).to_csv(os.path.join(working_directory, 
                                                                                          'cluster_representative.csv'), 
                                                                             index=False)
    
def get_expected_coverage(genome_length, reads_mapped, genome_totol_count):
    mean_mapping_length = genome_totol_count/reads_mapped
    
    N = genome_length/mean_mapping_length
    x = reads_mapped
    
    expected_M = N*(1-((1-1/N)**x))
    variance = N*((1-1/N)**x) + (N**2)*(1-1/N)*((1-2/N)**x)-(N**2)*((1-1/N)**(2*x))
    
    expected_coverage = expected_M/N
    try:
        std = math.sqrt(variance)
    except ValueError:
        std = 0
    return expected_coverage, std

def calculate_depth(assembly_id, taxa_records):
    genome_length = taxa_records[assembly_id]['genome_length']
    covbases = taxa_records[assembly_id]['covbases']
    genome_totol_count = taxa_records[assembly_id]['genome_totol_count']
    reads_mapped = taxa_records[assembly_id]['reads_mapped']
    
    if genome_totol_count == 0 or reads_mapped == 0:
        breadth_coverage = 0
        depth_coverage = 0
        expected_breadth_coverage = 0
    else:
        breadth_coverage = covbases/genome_length
        depth_coverage = genome_totol_count/covbases
        expected_breadth_coverage, std = get_expected_coverage(genome_length, reads_mapped, genome_totol_count)
    
    return breadth_coverage, depth_coverage, expected_breadth_coverage

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
parser.add_argument("--kingdom", type=str, help="A comma separated list of taxids of valid kingdoms. Default:[2,4751,2157,10239]", default='2,4751,2157,10239')
parser.add_argument("--include-mag", action='store_true', required=False, help="Include metagenomic assemble genomes. Default:[off]")
parser.set_defaults(include_mag=False)
parser.add_argument("--subspecies", action='store_true', required=False, help="Verify taxonomic classification at subspecies rank. Default:[off]")
parser.set_defaults(subspecies=False)
parser.add_argument("--accession", action='store_true', required=False, help="Take accession ids as taxids. Does not work with min-abundance. Default:[off]")
parser.set_defaults(accession=False)

# mob_output_paths = ['/home/Users/ns58/Mob-experiments/Mob-outputs/Zymo-Even',
#                     '/home/Users/ns58/Mob-experiments/Mob-outputs/Zymo-Even-0.25',
#                     '/home/Users/ns58/Mob-experiments/Mob-outputs/Zymo-Even-0.50',
#                     '/home/Users/ns58/Mob-experiments/Mob-outputs/Zymo-Even-0.75']

# input_fastq_paths =['/home/Users/ns58/Mob-experiments/Data-grouped/Zymo-EVEN-p0.01',
#                     '/home/Users/ns58/Mob-experiments/Data-grouped/Zymo-EVEN-p0.25',
#                     '/home/Users/ns58/Mob-experiments/Data-grouped/Zymo-EVEN-p0.50',
#                     '/home/Users/ns58/Mob-experiments/Data-grouped/Zymo-EVEN-p0.75']

# magnet_output = '/home/Users/yl181/memu/magnet_v222_ZYMO_EVEN'
# if not os.path.exists(magnet_output):
#     os.mkdir(magnet_output)
    
# cmd_records = []

# for idx, mob_output_path in enumerate(mob_output_paths):
#     for f in os.listdir(mob_output_path):
#         if f.startswith('Zymo-Even-') and 'relative_abundance.tsv' in f:
#             sample_id = '-'.join(f.split('-')[0:3])
#             input_fastq = os.path.join(input_fastq_paths[idx], f'{sample_id}.fastq')
#             if not os.path.exists(input_fastq):
#                 print("input_fastq missing:", mob_output_path, sample_id)
#             classification_output = os.path.join(mob_output_path, f)
#             if not os.path.exists(classification_output):
#                 print("classification_output missing:", input_fastq_paths[idx], sample_id)
                
#             magnet_output_prefix = input_fastq_paths[idx].split('/')[-1]
#             if not os.path.exists(os.path.join(magnet_output, magnet_output_prefix)):
#                 os.mkdir(os.path.join(magnet_output, magnet_output_prefix))
#             magnet_sample_output = os.path.join(magnet_output, magnet_output_prefix, sample_id)
            
#             # print(input_fastq)
#             # print(classification_output)
#             # print(magnet_sample_output)
            
#             cmd_args = ['-c', classification_output,
#             '-i', input_fastq,
#             '-m', 'ont',
#             '-a', str(12),
#             '--min-abundance', str(0.000),
#             '-o', magnet_sample_output,
#             '--threads', '40']
            
#             print(" ".join(cmd_args))
#             cmd_records.append(cmd_args)

# mob_output_paths = ['/home/Users/ns58/Mob-experiments/Mob-outputs/Cheese']
# input_fastq_paths = ['/home/Users/ns58/Mob-experiments/Data-grouped/Cheese']

# magnet_output = '/home/Users/yl181/memu/magnet_v222_cheese'
# if not os.path.exists(magnet_output):
#     os.mkdir(magnet_output)
    
# cmd_records = []

# for idx, mob_output_path in enumerate(mob_output_paths):
#     for f in os.listdir(mob_output_path):
#         if f.startswith('Cheese-') and 'relative_abundance.tsv' in f:
#             sample_id = '-'.join(f.split('-')[0:2])
#             input_fastq = os.path.join(input_fastq_paths[idx], f'{sample_id}.fastq')
#             if not os.path.exists(input_fastq):
#                 print("input_fastq missing:", mob_output_path, sample_id)
#             classification_output = os.path.join(mob_output_path, f)
#             if not os.path.exists(classification_output):
#                 print("classification_output missing:", input_fastq_paths[idx], sample_id)
                
#             magnet_output_prefix = input_fastq_paths[idx].split('/')[-1]
#             if not os.path.exists(os.path.join(magnet_output, magnet_output_prefix)):
#                 os.mkdir(os.path.join(magnet_output, magnet_output_prefix))
#             magnet_sample_output = os.path.join(magnet_output, magnet_output_prefix, sample_id)
            
#             # print(input_fastq)
#             # print(classification_output)
#             # print(magnet_sample_output)
            
#             cmd_args = ['-c', classification_output,
#             '-i', input_fastq,
#             '-m', 'ont',
#             '-a', str(12),
#             '--min-abundance', str(0.000),
#             '-o', magnet_sample_output,
#             '--threads', '40']
            
#             print(" ".join(cmd_args))
#             cmd_records.append(cmd_args)

# mob_output_paths = ['/home/Users/ns58/Mob-experiments/Mob-outputs/Zymo-Log-0.75']
# input_fastq_paths = ['/home/Users/ns58/Mob-experiments/Data-grouped/Zymo-Log-p0.75']

# magnet_output = '/home/Users/yl181/memu/magnet_v222_ZYMO_Log'
# if not os.path.exists(magnet_output):
#     os.mkdir(magnet_output)
    
# cmd_records = []

# for idx, mob_output_path in enumerate(mob_output_paths):
#     for f in os.listdir(mob_output_path):
#         if f.startswith('Zymo-Log-') and 'relative_abundance.tsv' in f:
#             sample_id = '-'.join(f.split('-')[0:3])
#             input_fastq = os.path.join(input_fastq_paths[idx], f'{sample_id}.fastq')
#             if not os.path.exists(input_fastq):
#                 print("input_fastq missing:", mob_output_path, sample_id)
#             classification_output = os.path.join(mob_output_path, f)
#             if not os.path.exists(classification_output):
#                 print("classification_output missing:", input_fastq_paths[idx], sample_id)
                
#             magnet_output_prefix = input_fastq_paths[idx].split('/')[-1]
#             if not os.path.exists(os.path.join(magnet_output, magnet_output_prefix)):
#                 os.mkdir(os.path.join(magnet_output, magnet_output_prefix))
#             magnet_sample_output = os.path.join(magnet_output, magnet_output_prefix, sample_id)
            
#             # print(input_fastq)
#             # print(classification_output)
#             # print(magnet_sample_output)
            
#             cmd_args = ['-c', classification_output,
#             '-i', input_fastq,
#             '-m', 'ont',
#             '-a', str(12),
#             '--min-abundance', str(0.000),
#             '-o', magnet_sample_output,
#             '--threads', '40']
            
#             print(" ".join(cmd_args))
#             cmd_records.append(cmd_args)

# mob_output_paths = ['/home/Users/ns58/Mob-experiments/Mob-outputs/Sim-MAG-Even-s20']
# input_fastq_paths = ['/home/Users/ns58/Mob-experiments/Data-grouped/Sim-MAG-Even-s20']

# magnet_output = '/home/Users/yl181/memu/magnet_v222_simulations_mags'
# if not os.path.exists(magnet_output):
#     os.mkdir(magnet_output)
    
# cmd_records = []

# for idx, mob_output_path in enumerate(mob_output_paths):
#     for f in os.listdir(mob_output_path):
#         if 'relative_abundance.tsv' in f:
#             sample_id = '-'.join(f.split('-')[0:3])
#             input_fastq = os.path.join(input_fastq_paths[idx], f'{sample_id}.fastq')
#             if not os.path.exists(input_fastq):
#                 print("input_fastq missing:", mob_output_path, sample_id)
#             classification_output = os.path.join(mob_output_path, f)
#             if not os.path.exists(classification_output):
#                 print("classification_output missing:", input_fastq_paths[idx], sample_id)
                
#             magnet_output_prefix = input_fastq_paths[idx].split('/')[-1]
#             if not os.path.exists(os.path.join(magnet_output, magnet_output_prefix)):
#                 os.mkdir(os.path.join(magnet_output, magnet_output_prefix))
#             magnet_sample_output = os.path.join(magnet_output, magnet_output_prefix, sample_id)
            
#             cmd_args = ['-c', classification_output,
#             '-i', input_fastq,
#             '-m', 'ont',
#             '-a', str(12),
#             '--min-abundance', str(0.000),
#             '-o', magnet_sample_output,
#             '--threads', '40']
            
#             #print(sample_id)
#             print(" ".join(cmd_args))
#             cmd_records.append(cmd_args)
            
# mob_output_paths = ['/home/Users/ns58/Mob-experiments/Mob-outputs/Zymo-TM']
# input_fastq_paths = ['/home/Users/ns58/Mob-experiments/Data-grouped/Zymo-TM']

# magnet_output = '/home/Users/yl181/memu/magnet_v222_ZYMO_TM'
# if not os.path.exists(magnet_output):
#     os.mkdir(magnet_output)
    
# cmd_records = []

# for idx, mob_output_path in enumerate(mob_output_paths):
#     for f in os.listdir(mob_output_path):
#         if 'relative_abundance.tsv' in f and not 'nof' in f:
#             sample_id = '-'.join(f.split('-')[0:3])
#             input_fastq = os.path.join(input_fastq_paths[idx], f'{sample_id}.fastq')
#             if not os.path.exists(input_fastq):
#                 print("input_fastq missing:", mob_output_path, sample_id)
#             classification_output = os.path.join(mob_output_path, f)
#             if not os.path.exists(classification_output):
#                 print("classification_output missing:", input_fastq_paths[idx], sample_id)
                
#             magnet_output_prefix = input_fastq_paths[idx].split('/')[-1]
#             if not os.path.exists(os.path.join(magnet_output, magnet_output_prefix)):
#                 os.mkdir(os.path.join(magnet_output, magnet_output_prefix))
#             magnet_sample_output = os.path.join(magnet_output, magnet_output_prefix, sample_id)
            
#             cmd_args = ['-c', classification_output,
#             '-i', input_fastq,
#             '-m', 'ont',
#             '-a', str(12),
#             '--min-abundance', str(0.000),
#             '-o', magnet_sample_output,
#             '--threads', '40']
            
#             #print(sample_id)
#             print(" ".join(cmd_args))
#             cmd_records.append(cmd_args)
            
# for idx, mob_output_path in enumerate(mob_output_paths):
#     for f in os.listdir(mob_output_path):
#         if 'relative_abundance.tsv' in f and 'nof' in f:
#             sample_id = '-'.join(f.split('-')[1:4])
#             input_fastq = os.path.join(input_fastq_paths[idx], f'{sample_id}.fastq')
#             if not os.path.exists(input_fastq):
#                 print("input_fastq missing:", mob_output_path, sample_id)
#             classification_output = os.path.join(mob_output_path, f)
#             if not os.path.exists(classification_output):
#                 print("classification_output missing:", input_fastq_paths[idx], sample_id)
                
#             magnet_output_prefix = '-'.join(f.split('-')[0:3])
#             if not os.path.exists(os.path.join(magnet_output, magnet_output_prefix)):
#                 os.mkdir(os.path.join(magnet_output, magnet_output_prefix))
#             magnet_sample_output = os.path.join(magnet_output, magnet_output_prefix, sample_id)
            
#             cmd_args = ['-c', classification_output,
#             '-i', input_fastq,
#             '-m', 'ont',
#             '-a', str(12),
#             '--min-abundance', str(0.000),
#             '-o', magnet_sample_output,
#             '--threads', '40']
            
#             #print(sample_id)
#             print(" ".join(cmd_args))
#             cmd_records.append(cmd_args)
            
    
# mob_output_paths = ['/home/Users/ns58/Mob-experiments/Mob-outputs/SRR17687125']
# input_fastq_paths = ['/home/Users/ns58/Mob-experiments/Data/Gut']

# magnet_output = '/home/Users/yl181/memu/magnet_v222_SRR17687125'
# if not os.path.exists(magnet_output):
#     os.mkdir(magnet_output)
    
# cmd_records = []

# for idx, mob_output_path in enumerate(mob_output_paths):
#     for f in os.listdir(mob_output_path):
#         if 'relative_abundance.tsv' in f and not 'nof' in f:
#             sample_id = f.split('-')[0]
#             input_fastq = os.path.join(input_fastq_paths[idx], f'{sample_id}.fastq')
#             if not os.path.exists(input_fastq):
#                 print("input_fastq missing:", mob_output_path, sample_id)
#             classification_output = os.path.join(mob_output_path, f)
#             if not os.path.exists(classification_output):
#                 print("classification_output missing:", input_fastq_paths[idx], sample_id)
                
#             magnet_output_prefix = input_fastq_paths[idx].split('/')[-1]
#             if not os.path.exists(os.path.join(magnet_output, magnet_output_prefix)):
#                 os.mkdir(os.path.join(magnet_output, magnet_output_prefix))
#             magnet_sample_output = os.path.join(magnet_output, magnet_output_prefix, sample_id)
            
#             cmd_args = ['-c', classification_output,
#             '-i', input_fastq,
#             '-m', 'ont',
#             '-a', str(12),
#             '--min-abundance', str(0.000),
#             '-o', magnet_sample_output,
#             '--threads', '40']
            
#             #print(sample_id)
#             print(" ".join(cmd_args))
#             cmd_records.append(cmd_args)
            
# for idx, mob_output_path in enumerate(mob_output_paths):
#     for f in os.listdir(mob_output_path):
#         if 'relative_abundance.tsv' in f and 'nof' in f:
#             sample_id = f.split('-')[1]
#             input_fastq = os.path.join(input_fastq_paths[idx], f'{sample_id}.fastq')
#             if not os.path.exists(input_fastq):
#                 print("input_fastq missing:", mob_output_path, sample_id)
#             classification_output = os.path.join(mob_output_path, f)
#             if not os.path.exists(classification_output):
#                 print("classification_output missing:", input_fastq_paths[idx], sample_id)
                
#             magnet_output_prefix = '-'.join(f.split('-')[0:2])
#             if not os.path.exists(os.path.join(magnet_output, magnet_output_prefix)):
#                 os.mkdir(os.path.join(magnet_output, magnet_output_prefix))
#             magnet_sample_output = os.path.join(magnet_output, magnet_output_prefix, sample_id)
            
#             cmd_args = ['-c', classification_output,
#             '-i', input_fastq,
#             '-m', 'ont',
#             '-a', str(12),
#             '--min-abundance', str(0.000),
#             '-o', magnet_sample_output,
#             '--threads', '40']
            
#             #print(sample_id)
#             print(" ".join(cmd_args))
#             cmd_records.append(cmd_args)
            
# mob_output_paths = ['/home/Users/ns58/Mob-experiments/Kraken2-outputs/Zymo-Log-0.10']
# input_fastq_paths = ['/home/Users/ns58/Mob-experiments/Data-grouped/Zymo-Log-p0.10']

# magnet_output = '/home/Users/yl181/memu/magnet_v222_ZYMO_Log'
# if not os.path.exists(magnet_output):
#     os.mkdir(magnet_output)
    
# cmd_records = []

# for idx, mob_output_path in enumerate(mob_output_paths):
#     for f in os.listdir(mob_output_path):
#         if f.startswith('Zymo-Log-') and 'report' in f:
#             sample_id = f.split('.')[0]
#             input_fastq = os.path.join(input_fastq_paths[idx], f'{sample_id}.fastq')
#             if not os.path.exists(input_fastq):
#                 print("input_fastq missing:", mob_output_path, sample_id)
#             classification_output = os.path.join(mob_output_path, f)
#             if not os.path.exists(classification_output):
#                 print("classification_output missing:", input_fastq_paths[idx], sample_id)
                
#             magnet_output_prefix = input_fastq_paths[idx].split('/')[-1]+'-kraken'
#             if not os.path.exists(os.path.join(magnet_output, magnet_output_prefix)):
#                 os.mkdir(os.path.join(magnet_output, magnet_output_prefix))
#             magnet_sample_output = os.path.join(magnet_output, magnet_output_prefix, sample_id)
            
#             # print(input_fastq)
#             # print(classification_output)
#             # print(magnet_sample_output)
            
#             cmd_args = ['-c', classification_output,
#             '-i', input_fastq,
#             '-m', 'ont',
#             '-t', str(4),
#             '-a', str(0),
#             '--min-abundance', str(0.01),
#             '-o', magnet_sample_output,
#             '--threads', '40']
            
#             print(" ".join(cmd_args))
#             cmd_records.append(cmd_args)

cmd_records = []

magnet_output = '/home/Users/yl181/memu/magnet_v222_MetaMapsSim'
if not os.path.exists(magnet_output):
    os.mkdir(magnet_output)

classification_output = '/home/Users/ns58/Mob-experiments/Lemur-outputs/MetaMaps-sim/nof-i100-relative_abundance.tsv'
if not os.path.exists(classification_output):
    print("classification_output missing:", classification_output)
input_fastq = '/home/Users/ns58/Mob-experiments/Data/Metamaps-sim/simulations_i100_specifiedFrequencies/0/reads.fastq'
if not os.path.exists(input_fastq):
    print("input_fastq missing:", input_fastq)
    
magnet_output_prefix = 'nof-i100'
if not os.path.exists(os.path.join(magnet_output, magnet_output_prefix)):
    os.mkdir(os.path.join(magnet_output, magnet_output_prefix))
magnet_sample_output = os.path.join(magnet_output, magnet_output_prefix, magnet_output_prefix)

cmd_args = ['-c', classification_output,
'-i', input_fastq,
'-m', 'ont',
'-a', str(12),
'--min-abundance', str(0.000),
'-o', magnet_sample_output,
'--threads', '40']

print(" ".join(cmd_args))
cmd_records.append(cmd_args)

for cmd_args in cmd_records:
    run_magnet(cmd_args)