#!/bin/bash -ue
minimap2 -t 1 --secondary=no -ax map-ont -y polished_assembly_10.fasta intermediate.fastq > results/intermediate.sam

# Clean up intermediate files
rm intermediate.fastq
