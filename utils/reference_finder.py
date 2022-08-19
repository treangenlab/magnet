import os
import subprocess

from collections import defaultdict
from Bio import Entrez
import time

import pandas as pd
import json

def run_datasets_summary(taxid, flags, assembly_level='complete_genome'):
    args = ['datasets', 'summary', 
            'genome',
            'taxon', str(taxid),
            '--assembly-level', assembly_level,
            '--limit', str(1)]
    
    args = args + flags
    
    grepOut  = subprocess.run(args,
                              universal_newlines=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.DEVNULL)
    
    return json.loads(grepOut.stdout.strip())

def run_datasets_download(taxid, assembly_accession, working_dir):
    grepOut  = subprocess.run(['datasets', 'download', 
                               'genome',
                               'accession', str(assembly_accession),
                               '--exclude-genomic-cds',
                               '--exclude-gff3',
                               '--exclude-protein',
                               '--exclude-rna',
                               '--filename', f"{working_dir}/{taxid}.zip",
                               '--no-progressbar'],
                              check=True, stderr=subprocess.DEVNULL)
    
    return grepOut.returncode

def download_reference_genome(taxid, working_dir):
    # check reference and representative genomes first
    #source = 'RefSeq'
    assembly_accession = None
    res =None
    representative = True # assembly_category
    assembly_level = 'complete'
    
    res = run_datasets_summary(taxid, ['--reference', '--assembly-source', 'refseq'])
    
    # check complete assembly in RefSeq
    if res['total_count'] == 0:
        representative = False
        res = run_datasets_summary(taxid, ['--assembly-source', 'refseq'])

    # check complete assembly
    if res['total_count'] == 0:
        res = run_datasets_summary(taxid, [])
        #source = 'Assembly'
    
    # check chromosome assembly
    if res['total_count'] == 0:
        assembly_level = 'chromosome'
        res = run_datasets_summary(taxid, [], assembly_level)
        
    # check chromosome assembly
    if res['total_count'] == 0:
        assembly_level='contig'
        res = run_datasets_summary(taxid, [], assembly_level)

    # check chromosome assembly
    if res['total_count'] == 0:
        assembly_level='scaffold'
        res = run_datasets_summary(taxid, [], assembly_level)
        
    if res['total_count'] == 1:
        assembly_accession = res['assemblies'][0]['assembly']['assembly_accession']
        try:
            source = res['assemblies'][0]['assembly']['annotation_metadata']['source']
        except KeyError:
            source = None
        assembly_level_ret = res['assemblies'][0]['assembly']['assembly_level']
        organism = res['assemblies'][0]['assembly']['org']['sci_name']
        strain = res['assemblies'][0]['assembly']['org']['strain']
        
        print(str(taxid).ljust(10), '\t', assembly_accession, '\t', organism, '\t', strain, '\t', assembly_level_ret)
        
        return_code = run_datasets_download(taxid, assembly_accession, working_dir)
        
        if not return_code:
            download_completed = True
        else:
            download_completed = False
            
    else:
        print(str(taxid).ljust(10), '\t', 'Genome Not Found.')
        assembly_accession = None
        source = None
        representative = None
        assembly_level_ret = None
        organism = None
        strain = None
        download_completed = False
        
    return taxid, assembly_accession, source, representative, assembly_level_ret, organism, strain, download_completed

def unpack(working_dir, output_dir):
    subprocess.run(['unzip', '-o', '-q', 
					f'{working_dir}/*.zip', 
					'-d', output_dir],
				   check=True)

def cat_reference_genome(reference_metadata, output_dir, reference_genome_path='reference_genomes'):
    if not os.path.exists(reference_genome_path):
        os.mkdir(reference_genome_path)
        
    downloaded_assemblies = reference_metadata[reference_metadata['Downloaded']]
    
    for assembly_id in downloaded_assemblies['Assembly Accession ID']:
        genome_path = os.path.join(output_dir, 'ncbi_dataset', 'data', assembly_id)
        output_fasta = os.path.join(reference_genome_path, f'{assembly_id}.fasta')

        subprocess.Popen(f'cat {genome_path}/*.fna', shell=True, stdout=open(os.path.join(output_fasta), "w"))
