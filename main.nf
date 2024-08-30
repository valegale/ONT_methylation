#!/usr/bin/env nextflow

include { bam2fastq; minimap2; sort; index } from './modules/map_index_bam.nf'
include { modkit_pileup } from './modules/modkit.nf'

workflow {
    bam_file = file('../v5_evaluation/barcode10/calling_barcode10.bam')
    reference = file('../v5_evaluation/barcode10/polished_assembly_10.fasta')
    output_dir = file('results')

    // Preprocessing for modkit
    fastq_file = bam2fastq(bam_file, output_dir)
    sam_file = minimap2(fastq_file, reference, output_dir)
    sorted_bam = sort(sam_file, output_dir)
    index(sorted_bam)

    // run modkit pileup 
    modkit_pileup(sorted_bam, output_dir)
}