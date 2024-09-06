#!/usr/bin/env nextflow

include { bam2fastq; minimap2; index } from './modules/map_index_bam.nf'
include { modkit_pileup; modkit_find_motifs } from './modules/modkit.nf'

workflow {
    bam_file = file('../v5_evaluation/barcode10/calling_barcode10.bam')
    reference = file('../v5_evaluation/barcode10/polished_assembly_10.fasta')

    // Preprocessing for modkit
    fastq_file = bam2fastq(bam_file)
    sam_file = minimap2(fastq_file, reference)
    index_bam = index(sam_file)

    // run modkit pileup 
    modkit_bed = modkit_pileup(sam_file, index_bam)

    // find motifs with modkit
    //modkit_find_motifs(modkit_bed, reference)
}