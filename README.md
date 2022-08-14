# Reference Finder

This pipeline handles reference genome selection, download using NCBI Datasets, read alignment using Minimap2, and breadth and depth of genome coverage calculations.

## Dependencies

- Python 3.9
- NCBI Datasets v13.30.2
- Minimap2 v2.24-r1122
- Samtools v1.15.1
- Biopython
- Pandas
- Ete3

## Getting started

```
usage: reference_inference.py [-h] -i FASTA -o OUTPUT -w WORKING -d DATABASES [-m MIN_FRAC] [-c MIN_COVERAGE_SCORE]
                              [--online] [-t THREADS]

Repeated Read Alignment Module

optional arguments:
  -h, --help            show this help message and exit
  -i FASTA, --fasta FASTA
                        Sequences in fasta format.
  -o OUTPUT, --output OUTPUT
                        Output directory.
  -w WORKING, --working WORKING
                        Working directory.
  -d DATABASES, --databases DATABASES
                        SeqScreenDB
  -m MIN_FRAC, --min-frac MIN_FRAC
                        minimum fraction of assigned reads for a species to be included in the first alignment
                        process. [0.002]
  -c MIN_COVERAGE_SCORE, --min-coverage-score MIN_COVERAGE_SCORE
                        minimum coverage score for a species to be included in the second alignment process. [0.7]
  --online              Use online mode for searching reference genome. Requires internet access.
  -t THREADS, --threads THREADS
                        Number of threads. [1]
```

## Example Output

The following table is an example of the output csv file.

| Taxonomy ID | Assembly Accession ID | Source Database | Is Representative | Assembly Level | Organism of Assembly | Downloaded | Species | Breadth Coverage | Expected Coverage | Coverage Score | Depth Coverage | BC2  | EC2  | CS2  | DC2  | Consensus ANI |
|-------------|-----------------------|-----------------|-------------------|----------------|----------------------|------------|---------|------------------|-------------------|----------------|----------------|------|------|------|------|---------------|
| 1747    | GCF_000376705.1 | NCBI RefSeq                                   | TRUE  | Complete Genome | Cutibacterium acnes HL096PA1                  | TRUE | Cutibacterium acnes     | 0.41 | 0.47 | 0.87 | 1.54 | 0.41 | 0.4  | 1.03 | 1.24 | 0.79 | 
| 1280    | GCF_000013425.1 | University of Oklahoma Health Sciences Center | TRUE  | Complete Genome | Staphylococcus aureus subsp. aureus NCTC 8325 | TRUE | Staphylococcus aureus   | 0.87 | 0.92 | 0.95 | 2.87 | 0.86 | 0.87 | 0.99 | 2.35 | 0.97 | 
| 210     | GCF_017821535.1 | NCBI RefSeq                                   | TRUE  | Complete Genome | Helicobacter pylori                           | TRUE | Helicobacter pylori     | 0.63 | 0.72 | 0.88 | 2.01 | 0.63 | 0.65 | 0.97 | 1.66 | 0.85 | 
| 1396    | GCF_002220285.1 | NCBI RefSeq                                   | TRUE  | Complete Genome | Bacillus cereus                               | TRUE | Bacillus cereus         | 0.74 | 0.96 | 0.77 | 4.26 | 0.02 | 0.05 | 0.46 | 2.25 | 0.84 | 
| 1890302 | GCF_008807735.1 | NCBI RefSeq                                   | TRUE  | Complete Genome | Bacillus wiedmannii                           | TRUE | Bacillus wiedmannii     | 0.77 | 0.96 | 0.79 | 4.32 | 0.04 | 0.08 | 0.45 | 2.31 | 0.82 | 
| 1423    | GCF_000009045.1 | BSNR                                          | TRUE  | Complete Genome | Bacillus subtilis subsp. subtilis str. 168    | TRUE | Bacillus subtilis       | 0.03 | 0.23 | 0.13 | 8.75 | 0    | 0    | 0    | 0    | 0    | 
| 294     | GCF_000730425.1 | NCBI RefSeq                                   | FALSE | Complete Genome | Pseudomonas fluorescens                       | TRUE | Pseudomonas fluorescens | 0.21 | 0.36 | 0.57 | 2.17 | 0    | 0    | 0    | 0    | 0    | 

The table contains the following feature of about the reference genome that are being selected for alignment:

- Taxonomy ID: The taxonomy ID of the species that users queries
- Assembly Accession ID: The assembly accession ID of the genome that is being selected as reference genome for the alignment process
- Source Database: source of the assembly, part of the metadata from NCBI Database
- Is Representative: Whether or not the assembly is labeled as representative genome in NCBI database
- Assembly Level: the assembly level of the assembly, can be one of the following: complete genome, chromosome, contig, scaffold
- Organism of Assembly: The description of the assembly, may contain strain level information
- Species: taxonomy name of the species
