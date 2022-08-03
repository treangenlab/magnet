from concurrent.futures import thread
import os
import subprocess
import pandas as pd

def sort_samfile(assembly_id, output_dir, num_cores):
    '''converting and sorting alignment files'''
    bam_files = os.path.join(output_dir, "bam_files")
    sam_files = os.path.join(output_dir, "sam_files")
    
    if not os.path.exists(bam_files):
        os.mkdir(bam_files)

    # covert sam file to binary bam file
    subprocess.run([
        "samtools",
        "view",
        "-@", str(num_cores),
        "-bS", os.path.join(sam_files, f"{assembly_id}.sam"),
        "-o", os.path.join(bam_files, f"{assembly_id}.bam")],
                    check=True)

    # sort the bam file 
    subprocess.run([
        "samtools",
        "sort",
        "-@", str(num_cores),
        "-o", os.path.join(bam_files, f"{assembly_id}.sorted.bam"),
        "-O", "BAM", os.path.join(bam_files, f"{assembly_id}.bam")],
                    check=True)

    # indexing the sorted bam file 
    subprocess.run([
        "samtools",
        "index",
        os.path.join(bam_files, f"{assembly_id}.sorted.bam")],
                    check=True)

def run_minimap2(input_fastq, reference_file, assembly_id, output_dir, threads=20):

    sam_files = os.path.join(output_dir, "sam_files")

    if not os.path.exists(sam_files):
        os.mkdir(sam_files)

    subprocess.run(["minimap2", 
                    "-ax", "map-ont", 
                    reference_file, 
                    input_fastq, 
                    "--sam-hit-only", 
                    "-o", os.path.join(sam_files, f"{assembly_id}.sam"),
                    "-t", str(threads)],
                    check=True)

def calculate_depth(assembly_id, output_dir):
    depth_files = os.path.join(output_dir, "depth_files")
    bam_files = os.path.join(output_dir, "bam_files")

    if not os.path.exists(depth_files):
        os.mkdir(depth_files)

    depth_file = os.path.join(depth_files, f"{assembly_id}.depth")
                
    subprocess.run([
        "samtools",
        "depth",
        os.path.join(bam_files, f"{assembly_id}.sorted.bam")],
                    check=True,
                    stdout=open(os.path.join(depth_file), "w"))
        

reference_metadata = pd.read_csv('reference_metadata.csv')
reference_genome_path = 'reference_genomes'
input_fastq = '/home/Users/yl181/seqscreen_nano/ZymoBIOMICS.STD2.Log.ont.raw_sequences/ERR3152366.downsampled.fastq'

downloaded_assemblies = reference_metadata[reference_metadata['Downloaded']]

output_dir = 'alignment'
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

num_cores = 20

for assembly_id in downloaded_assemblies['Assembly Accession ID']:
    reference_fasta = os.path.join(reference_genome_path, f'{assembly_id}.fasta')
    #run_minimap2(input_fastq, reference_fasta, assembly_id, output_dir, threads=num_cores)
    #sort_samfile(assembly_id, output_dir, num_cores)
    calculate_depth(assembly_id, output_dir)
    