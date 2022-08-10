# Reference Finder

This pipeline handles reference genome selection, download using NCBI Datasets, read alignment using Minimap2, and breadth and depth of genome coverage calculations.

## Getting started

```
usage: pipeline.py [-h] -i FASTA -o OUTPUT [-w WORKING] [-m MIN_FRAC] [-c MIN_COVERAGE_SCORE] [-t THREADS]

Repeated Read Alignment Module

optional arguments:
  -h, --help            show this help message and exit
  -i FASTA, --fasta FASTA
                        Sequences in fasta format.
  -o OUTPUT, --output OUTPUT
                        seqscreen output directory.
  -w WORKING, --working WORKING
                        working directory. [working]
  -m MIN_FRAC, --min-frac MIN_FRAC
                        minimum fraction of assigned reads for a species to be included in the first alignment
                        process. [0.002]
  -c MIN_COVERAGE_SCORE, --min-coverage-score MIN_COVERAGE_SCORE
                        minimum coverage score for a species to be included in the second alignment process. [0.7]
  -t THREADS, --threads THREADS
                        number of threads. [1]
```

## Dependencies

- Python 3.9
- NCBI Datasets v13.30.2
- Minimap2 v2.17-r941
- Samtools v1.7
- Biopython
- Pandas
- Ete3
