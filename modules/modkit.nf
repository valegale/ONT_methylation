process modkit_pileup {
    // execute the modkit pileup command
    publishDir  params.outdir, mode:'copy'

    cpus 8

    input:
    tuple val(sample_id), path(mapped_bam), path(index_bam)

    output:
    tuple val(sample_id), path("${sample_id}_modkit_pileup_output.bed")

    script:
    """
    modkit pileup -t $task.cpus  ${mapped_bam} ${sample_id}_modkit_pileup_output.bed --filter-threshold ${params.filter_threshold_modkit}
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
    modkit find-motifs -t 12 --in-bedmethyl ${bed_file} --ref ${reference} -o outout.tsv
    """
}

