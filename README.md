# Reference Finder

This pipeline handles reference genome selection, download using NCBI Datasets, read alignment using Minimap2, and breadth and depth of genome coverage calculations.

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

## Dependencies

- Python 3.9
- NCBI Datasets v13.30.2
- Minimap2 v2.24-r1122
- Samtools v1.15.1
- Biopython
- Pandas
- Ete3



