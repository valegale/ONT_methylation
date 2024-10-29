#!/usr/bin/env nextflow

include { bam2fastq; zipfastq; minimap2; index } from './modules/map_index_bam.nf'
include { modkit_pileup; modkit_pileup_bedgraphs; modkit_find_motifs; custom_bedgraphs} from './modules/modkit.nf'
include { compute_statistics } from './modules/statistics.nf'


    // Define the list of BAM files and reference filescd
    params.input_files = [

       [ file('../ONT_christine/results/BruAbortus.bam'), file('../ONT_christine/reference_genomes/BruAbortus_reference.fasta')]

       // [ file('../run_dorado_v5/barcode01/calling_barcode1.bam'), file('../run_dorado_v5/barcode01/polished_assembly_01.fasta')],
       // [ file('../run_dorado_v5/barcode02/calling_barcode02.bam'), file('../run_dorado_v5/barcode02/polished_assembly_02.fasta')],
       // [ file('../run_dorado_v5/barcode03/calling_barcode03.bam'), file('../run_dorado_v5/barcode03/polished_assembly_03.fasta')],
       // [ file('../run_dorado_v5/barcode04/calling_barcode4.bam'), file('../run_dorado_v5/barcode04/polished_assembly_04.fasta')],
       // [ file('../run_dorado_v5/barcode05/calling_barcode5.bam'), file('../run_dorado_v5/barcode05/polished_assembly_05.fasta')],
       // [ file('../run_dorado_v5/barcode06/calling_barcode6.bam'), file('../run_dorado_v5/barcode06/polished_assembly_06.fasta')],
       // [ file('../run_dorado_v5/barcode07/calling_barcode7.bam'), file('../run_dorado_v5/barcode07/polished_assembly_07.fasta')],
       // [ file('../run_dorado_v5/barcode08/calling_barcode08.bam'), file('../run_dorado_v5/barcode08/polished_assembly_08.fasta')],
       // [ file('../run_dorado_v5/barcode09/calling_barcode09.bam'), file('../run_dorado_v5/barcode09/polished_assembly_09.fasta')],
       // [ file('../run_dorado_v5/barcode10/calling_barcode10.bam'), file('../run_dorado_v5/barcode10/polished_assembly_10.fasta')],
       // [ file('../run_dorado_v5/barcode11/calling_barcode11.bam'), file('../run_dorado_v5/barcode11/polished_assembly_11.fasta')],
       // [ file('../run_dorado_v5/barcode12/calling_barcode12.bam'), file('../run_dorado_v5/barcode12/polished_assembly_12.fasta')]
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
    custom_bedgraphs(bed_file)
    compute_statistics(bed_file)

}


//  nextflow run main.nf -profile rki_slurm,rki_singularity -c rki_profile.config