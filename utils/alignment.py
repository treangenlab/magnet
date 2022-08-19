import os
import subprocess
import argparse
import sys
import pandas as pd

def sort_samfile(assembly_id, output_dir, num_cores):
    '''converting and sorting alignment files'''
    bam_files = os.path.join(output_dir, "bam_files")
    sam_files = os.path.join(output_dir, "sam_files")
    
    if not os.path.exists(bam_files):
        os.mkdir(bam_files)

    # covert sam file to binary bam file
    samtools_view_res = subprocess.Popen([
        "samtools",
        "view",
        "-@", str(num_cores),
        "-bS", os.path.join(sam_files, f"{assembly_id}.sam")],
        stdout=subprocess.PIPE)

    # sort the bam file 
    subprocess.run([
        "samtools",
        "sort",
        "-@", str(num_cores),
        "-o", os.path.join(bam_files, f"{assembly_id}.sorted.bam"),
        "-O", "BAM"],
        stdin=samtools_view_res.stdout,
        stderr=subprocess.DEVNULL,
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
                    check=True,
                    stderr=subprocess.DEVNULL)
        
def samtools_calculate_depth(assembly_id, output_dir):
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
def main(argv):
	'''main function'''
	parser = argparse.ArgumentParser(description="Read Alignment and Coverage Calculation.")
	parser.add_argument("-i", "--fastq", type=str, help="Sequences in fastq format", required=True)
	parser.add_argument("-o", "--output", type=str, help="Output directory", required=True)
	parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of threads. [1]")
	
	args = parser.parse_args()
	
	output_dir = args.output
	input_fastq = args.fastq
	threads = args.threads
	
	reference_metadata = pd.read_csv(os.path.join(output_dir, 'reference_metadata.csv'))
	reference_genome_path = os.path.join(output_dir, 'reference_genomes')

	downloaded_assemblies = reference_metadata[reference_metadata['Downloaded']]

	for assembly_id in downloaded_assemblies['Assembly Accession ID']:
		reference_fasta = os.path.join(reference_genome_path, f'{assembly_id}.fasta')
		run_minimap2(input_fastq, reference_fasta, assembly_id, output_dir, threads=threads)
		sort_samfile(assembly_id, output_dir, threads)
		calculate_depth(assembly_id, output_dir)

if __name__ == "__main__":
	main(sys.argv[1:])
