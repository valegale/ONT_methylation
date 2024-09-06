process modkit_pileup {
    // execute the modkit pileup command
    publishDir  params.outdir, mode:'copy'

    cpus 8

    input:
    path bam_file
    path bam_file_index

    output:
    path "modkit_pileup_output.bed" 

    script:
    """
    modkit pileup -t $task.cpus  ${bam_file} modkit_pileup_output.bed --filter-threshold ${params.filter_threshold_modkit}
    """
}

process modkit_find_motifs {
    // find motifs from the output of modkit pileup
    publishDir  params.outdir, mode:'copy'

    input:
    path bed_file
    path reference

    output:
    path "motifs.log" 

    script:
    """
    modkit find-motifs -t 12 --in-bedmethyl ${bed_file} --ref ${reference}  > motifs.log
    """
}

