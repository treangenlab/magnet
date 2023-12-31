import os
import subprocess
import time
import json

from collections import defaultdict
from Bio import Entrez
import pandas as pd

def run_datasets_summary(taxid, flags, assembly_level='complete'):
    args = ['datasets', 'summary', 
            'genome',
            'taxon', str(taxid),
            '--exclude-atypical',
            '--assembly-level', assembly_level,
            '--limit', '1']
    
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
                               '--include', 'genome',
                               '--filename', f"{working_dir}/{taxid}.zip",
                               '--no-progressbar'],
                              check=True, stderr=subprocess.DEVNULL)
    
    return grepOut.returncode

def cut_text(text, text_length):
    try:
        if len(text) > text_length-2:
            return text[0:text_length-2] + ".."
        else:
            return text.ljust(text_length)
    except TypeError:
        return "".ljust(text_length)
	
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
        
    # check contig assembly
    if res['total_count'] == 0:
        assembly_level='contig'
        res = run_datasets_summary(taxid, [], assembly_level)

    # check scaffold assembly
    if res['total_count'] == 0:
        assembly_level='scaffold'
        res = run_datasets_summary(taxid, [], assembly_level)
        
    if res['total_count'] > 0:
        assembly_accession = res['reports'][0]['accession']
        try:
            source = res['reports'][0]['source_database']
        except KeyError:
            source = None
        assembly_level_ret = res['reports'][0]['assembly_info']['assembly_level']
        organism = res['reports'][0]['organism']['organism_name']
        total_length = int(res['reports'][0]['assembly_stats']['total_sequence_length'])
        try:
            strain = res['reports'][0]['organism']['infraspecific_names']['strain']
        except KeyError:
            strain = None
		
        print(str(taxid).ljust(10), '\t', assembly_accession, '\t', cut_text(organism, 30), '\t', cut_text(strain, 10), '\t', assembly_level_ret)
        
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
        total_length = None
        download_completed = False
        
    return taxid, assembly_accession, source, representative, assembly_level_ret, organism, strain, total_length, download_completed

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

def prepare_reference_genomes(taxid_queries, output_directory, ncbi_taxa_db):
    working_dir = os.path.join(output_directory, 'ncbi_downloads')
    if not os.path.exists(working_dir):
        os.mkdir(working_dir)
        
    download_result = []
    for taxid in taxid_queries:
        download_result.append(download_reference_genome(taxid, working_dir))
        
    unpack(working_dir, output_directory)
    reference_metadata = pd.DataFrame(download_result,
                                      columns=['Taxonomy ID', 
                                               'Assembly Accession ID', 
                                               'Source Database', 
                                               'Is Representative', 
                                               'Assembly Level', 
                                               'Organism of Assembly',
                                               'Strain',
											   'Total Length',
                                               'Downloaded'])
    
    taxonomy_name = []
    for taxid in reference_metadata['Taxonomy ID']:
        taxonomy_name.append(ncbi_taxa_db.get_taxid_translator([taxid])[taxid])
    reference_metadata['Species'] = taxonomy_name
    
    reference_metadata.to_csv(os.path.join(output_directory, 'reference_metadata.csv'), index=False)
    cat_reference_genome(reference_metadata, output_directory, reference_genome_path=os.path.join(output_directory, 'reference_genomes'))
    
    return reference_metadata
