#!/usr/bin/env nextflow

include { bam2fastq; zipfastq; minimap2; index } from './modules/map_index_bam.nf'
include { modkit_pileup; modkit_find_motifs } from './modules/modkit.nf'


    // Define the list of BAM files and reference files
    params.input_files = [
        [ file('../ralstonia/barcode02/calling_barcode02.bam'), file('../ralstonia/assemblies/23_B33984_02_ONT_polished.fa') ],
        [ file('../ralstonia/barcode01/calling_barcode01.bam'), file('../ralstonia/assemblies/13_B3672_01_ONT_polished.fa') ]
        [ file('../ralstonia/barcode03/calling_barcode03.bam'), file('../ralstonia/assemblies/BK115556_03_ONT_polished.fa') ],
        [ file('../ralstonia/barcode03/calling_barcode03_duplex.bam'), file('../ralstonia/assemblies/BK115556_03_ONT_polished.fa') ],
        [ file('../ralstonia/barcode04/calling_barcode04.bam'), file('../ralstonia/assemblies/BB_XA27_0046_04_ONT_polished.fa') ],
        [ file('../ralstonia/barcode04/calling_barcode04_duplex.bam'), file('../ralstonia/assemblies/BB_XA27_0046_04_ONT_polished.fa') ],
        [ file('../ralstonia/barcode05/calling_barcode05.bam'), file('../ralstonia/assemblies/BB_XB07_0060_05_ONT_polished.fa') ],
        [ file('../ralstonia/barcode05/calling_barcode05_duplex.bam'), file('../ralstonia/assemblies/BB_XB07_0060_05_ONT_polished.fa') ]
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
    
    modkit_pileup(bam_and_index)
}
