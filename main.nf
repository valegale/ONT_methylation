#!/usr/bin/env nextflow

include { bam2fastq; zipfastq; minimap2; index } from './modules/map_index_bam.nf'
include { modkit_pileup; modkit_find_motifs } from './modules/modkit.nf'

workflow {
    bam_file = file('../ralstonia/barcode02/calling_barcode02.bam')
    reference = file('../ralstonia/assemblies/23_B33984_02_ONT_polished.fa')

    // Preprocessing for modkit
    fastq_file = bam2fastq(bam_file)
    zipfastq(fastq_file)
    sam_file = minimap2(fastq_file, reference)
    index_bam = index(sam_file)

    // run modkit pileup 
    modkit_bed = modkit_pileup(sam_file, index_bam)

    // find motifs with modkit
    //modkit_find_motifs(modkit_bed, reference)
}
