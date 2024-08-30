#!/bin/bash -ue
samtools view -b intermediate.sam |     samtools sort -@ 10 -o results/methylation_mapped.bam

# Clean up intermediate files
rm intermediate.sam
