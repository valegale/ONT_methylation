process bam2fastq {
    // converts the BAM file from Dorado to FASTQ format 
    input:
    tuple val(sample_id), path(bam_file), path(reference) 

    output:
    tuple val(sample_id), path("${sample_id}.fastq"), path(reference) 

    script:
    """
    samtools fastq ${bam_file} -T MM,ML > ${sample_id}.fastq

    """
}

process zipfastq {
    // zip fastq file
    publishDir  params.outdir, mode:'copy'
        
    input:
    tuple val(sample_id), path(fastq_file), path(reference)

    output:
    path "${sample_id}.fastq.gz"

    script:
    """
    gzip -c ${fastq_file} > ${sample_id}.fastq.gz
    """
}

process minimap2 {
    // align with minimap2 
    publishDir  params.outdir, mode:'copy'

    cpus 8
     
    input:
    tuple val(sample_id), path(fastq_file), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}_methylation_mapped.bam"), path(reference)

    script:
    """
    minimap2 -t $task.cpus --secondary=no -ax map-ont -y ${reference} ${fastq_file} | \
    samtools view -b | \
    samtools sort -@ $task.cpus -o ${sample_id}_methylation_mapped.bam
    """
}

process index {
    // index an aligned and sorted BAM
    publishDir  params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(mapped_bam), path(reference)

    output:
    tuple val(sample_id), path ("${mapped_bam}.bai")

    script:
    """
    samtools index ${mapped_bam}
    """
}
