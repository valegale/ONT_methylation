process modkit_pileup {
    // execute the modkit pileup command
    publishDir  params.outdir, mode:'copy'

    cpus 8

    input:
    tuple val(sample_id), path(mapped_bam), path(reference), path(index_bam)

    output:
    tuple val(sample_id), path("${sample_id}/modkit_pileup_output.bed"), path(reference)

    script:
    """
    modkit pileup -t $task.cpus  ${mapped_bam} ${sample_id}/modkit_pileup_output.bed --filter-threshold ${params.filter_threshold_modkit} 
    """
}

process modkit_pileup_bedgraphs {
    // execute the modkit pileup command
    publishDir  params.outdir, mode:'copy'

    cpus 8

    input:
    tuple val(sample_id), path(mapped_bam), path(reference), path(index_bam)

    output:
    path("${sample_id}/bedgraphs")

    script:
    """
    modkit pileup -t $task.cpus  ${mapped_bam} --bedgraph ${sample_id}/bedgraphs --filter-threshold ${params.filter_threshold_modkit} 
    """
}

process modkit_find_motifs {
    // find motifs from the output of modkit pileup
    publishDir  params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(bed_file), path(reference)

    output:
    path "${sample_id}/modkit_motifs.tsv" 

    script:
    """
    mkdir -p ${sample_id}
    modkit find-motifs -t 12 --in-bedmethyl ${bed_file} --ref ${reference} -o ${sample_id}/modkit_motifs.tsv
    """
}

