import os
import subprocess
import sys
import math
from collections import defaultdict
from Bio import Entrez
import time

import pandas as pd
import json
from ete3 import NCBITaxa
from Bio import SeqIO

def alignment_1_summary(downloaded_assemblies, output_directory):
    breadth_coverage_list = []
    depth_coverage_list = []
    expected_breadth_coverage_list = []
    coverage_score = []
    for assembly_id in downloaded_assemblies['Assembly Accession ID']:
        breadth_coverage, depth_coverage, expected_breadth_coverage = calculate_depth(assembly_id, output_directory, min_depth=1)
        breadth_coverage_list.append(breadth_coverage)
        depth_coverage_list.append(depth_coverage)
        expected_breadth_coverage_list.append(expected_breadth_coverage)
        if expected_breadth_coverage != 0:
            coverage_score.append(breadth_coverage/expected_breadth_coverage)
        else:
            coverage_score.append(0)
        
    downloaded_assemblies['Breadth Coverage'] = breadth_coverage_list
    downloaded_assemblies['Expected Coverage'] = expected_breadth_coverage_list
    downloaded_assemblies['Coverage Score'] = coverage_score
    downloaded_assemblies['Depth Coverage'] = depth_coverage_list
    
    downloaded_assemblies.to_csv(os.path.join(output_directory, 'alignment.csv'), index=False)
    
    return downloaded_assemblies

def calculate_depth(assembly_id, output_directory, min_depth=1):
    depth_file = os.path.join(output_directory, 'depth_files', f"{assembly_id}.depth")
    
    genome_length, genome_ids = parse_reference_fasta(assembly_id, output_directory)
    mapping_stats = get_mapping_stats(assembly_id, output_directory)
    reads_mapped = mapping_stats['reads mapped']

    genome_pos_count = 0
    genome_totol_count = 0
    with open(depth_file, "r") as depth:
        for line in depth.readlines():
            
            genome_id = line.split("\t")[0]
            pos = int(line.split("\t")[1])
            depth = int(line.strip().split("\t")[2])

            if depth >= min_depth and genome_id in genome_ids:
                genome_pos_count += 1
                genome_totol_count += depth

    
    if genome_totol_count == 0:
        breadth_coverage = 0
        depth_coverage = 0
        expected_breadth_coverage = 0
    else:
        breadth_coverage = genome_pos_count/genome_length
        depth_coverage = genome_totol_count/genome_pos_count
        expected_breadth_coverage, std = get_expected_coverage(genome_length, reads_mapped, genome_totol_count)
    
    return breadth_coverage, depth_coverage, expected_breadth_coverage

def parse_reference_fasta(assembly_id, output_directory):
    reference_fasta = os.path.join(output_directory, 'reference_genomes', f'{assembly_id}.fasta')
    
    with open(reference_fasta, "r") as handle:
        genome_length = 0
        genome_ids = []
        
        for record in SeqIO.parse(handle, "fasta"):      
            genome_length += len(record.seq)
            genome_ids.append(record.id)
    
    return genome_length, genome_ids

def get_mapping_stats(assembly_id, output_directory, genome_id=None):
    bam_file = os.path.join(output_directory, 'bam_files', f"{assembly_id}.sorted.bam")
    
    if genome_id is not None:
        stats_res = subprocess.Popen(['samtools', 'stats', bam_file, genome_id],
                               stdout=subprocess.PIPE)
    else:
        stats_res = subprocess.Popen(['samtools', 'stats', bam_file],
                               stdout=subprocess.PIPE)

    grep_res = subprocess.Popen(['grep', '^SN'],
                              stdin=stats_res.stdout,
                              stdout=subprocess.PIPE)

    mapping_res = subprocess.run(['cut', '-f', '2-'],
                                 check=True,
                                 universal_newlines=True,
                                 stdin = grep_res.stdout,
                                 stdout=subprocess.PIPE)

    mapping_stats = dict()
    for line in mapping_res.stdout.strip().split('\n'):
        attribute, value = line.split(':\t')
        mapping_stats[attribute] = float(value.split('\t')[0])
        
    return mapping_stats

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

def alignment_2_summary(downloaded_assemblies, output_directory):
    breadth_coverage_dict, depth_coverage_dict, expected_breadth_coverage_dict = \
    calculate_depth_merged(list(downloaded_assemblies['Assembly Accession ID']), output_directory, min_depth=1)
    
    breadth_coverage_list = []
    depth_coverage_list = []
    expected_breadth_coverage_list = []
    coverage_score = []
    for assembly_id in downloaded_assemblies['Assembly Accession ID']:
        breadth_coverage_list.append(breadth_coverage_dict[assembly_id])
        depth_coverage_list.append(depth_coverage_dict[assembly_id])
        expected_breadth_coverage_list.append(expected_breadth_coverage_dict[assembly_id])

        if expected_breadth_coverage_dict[assembly_id] != 0:
            coverage_score.append(breadth_coverage_dict[assembly_id]/expected_breadth_coverage_dict[assembly_id])
        else:
            coverage_score.append(0)
            
    downloaded_assemblies['BC2'] = breadth_coverage_list
    downloaded_assemblies['EC2'] = expected_breadth_coverage_list
    downloaded_assemblies['CS2'] = coverage_score
    downloaded_assemblies['DC2'] = depth_coverage_list
    
    downloaded_assemblies.sort_values(['CS2'], ascending=False).to_csv(os.path.join(output_directory, 'alignment.csv'), index=False)
    
    return downloaded_assemblies

def calculate_depth_merged(assembly_ids, output_directory, min_depth=1):
    depth_file = os.path.join(output_directory, 'depth_files', f"merged.depth")
    
    genome_id_pos_count = defaultdict(int)
    genome_id_totol_count = defaultdict(int)
    
    with open(depth_file, "r") as depth:
        for line in depth.readlines():            
            genome_id = line.split("\t")[0]
            pos = int(line.split("\t")[1])
            depth = int(line.strip().split("\t")[2])

            if depth >= min_depth:
                genome_id_pos_count[genome_id] += 1
                genome_id_totol_count[genome_id] += depth
    
    breadth_coverage_dict = defaultdict(float)
    depth_coverage_dict = defaultdict(float)
    reads_mapped_dict = defaultdict(int)
    expected_breadth_coverage_dict = defaultdict(float)
    
    for assembly_id in assembly_ids:
        genome_length, genome_ids = parse_reference_fasta(assembly_id, output_directory)
        
        genome_pos_count = 0
        genome_totol_count = 0
        for genome_id in genome_ids:
            genome_pos_count += genome_id_pos_count[genome_id]
            genome_totol_count += genome_id_totol_count[genome_id]
            if genome_id_pos_count[genome_id] > 0:
                mapping_stats = get_mapping_stats('merged', output_directory, genome_id)
                reads_mapped_dict[assembly_id] += mapping_stats['reads mapped']
            
        if genome_totol_count > 0:
            breadth_coverage_dict[assembly_id] = genome_pos_count/genome_length
            depth_coverage_dict[assembly_id] = genome_totol_count/genome_pos_count
            expected_breadth_coverage, std = get_expected_coverage(genome_length, reads_mapped_dict[assembly_id], genome_totol_count)
            expected_breadth_coverage_dict[assembly_id] = expected_breadth_coverage
            
    return breadth_coverage_dict, depth_coverage_dict, expected_breadth_coverage_dict

