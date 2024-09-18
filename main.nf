#!/usr/bin/env nextflow

include { bam2fastq; zipfastq; minimap2; index } from './modules/map_index_bam.nf'
include { modkit_pileup; modkit_pileup_bedgraphs; modkit_find_motifs } from './modules/modkit.nf'


    // Define the list of BAM files and reference files
    params.input_files = [
        [ file('../ralstonia/test_dataset/contig_4.bam'), file('../ralstonia/test_dataset/contig_4.fasta') ]
    ]

workflow {
    
    Channel
    .from( params.input_files )
    .map { bam_file, reference -> 
                sample_id = bam_file.simpleName 
                tuple(sample_id, bam_file, reference)  
            }
    .set { bam_ref_pairs }
    
    fastq_files = bam2fastq(bam_ref_pairs)
    zipfastq(fastq_files)
    mapped_bams = minimap2(fastq_files)
    index_bam = index(mapped_bams)

    // Join the channels by sample ID
    bam_and_index = mapped_bams
        .join(index_bam)
    
    bed_file = modkit_pileup(bam_and_index)
    modkit_pileup_bedgraphs(bam_and_index)
    modkit_find_motifs(bed_file)
}
