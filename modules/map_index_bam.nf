process bam2fastq {
    // converts the BAM file from Dorado to FASTQ format 
    input:
    path bam_file
    path output_dir

    output:
    path "${output_dir}/intermediate.fastq" 

    script:
    """
    samtools fastq ${bam_file} -T MM,ML > ${output_dir}/intermediate.fastq

    """
}

process minimap2 {
    // align with minimap2
    input:
    path fastq_file
    path reference
    path output_dir

    output:
    path "${output_dir}/intermediate.sam"

    script:
    """
    minimap2 -t $task.cpus --secondary=no -ax map-ont -y ${reference} ${fastq_file} > ${output_dir}/intermediate.sam

    # Clean up intermediate files
    rm ${fastq_file}
    """
}

process sort {
    // sort an aligned BAM
    input:
    path sam_file 
    path output_dir

    output:
    path "${output_dir}/methylation_mapped.bam" 

    script:
    """
    samtools view -b ${sam_file} | \
    samtools sort -@ 10 -o ${output_dir}/methylation_mapped.bam

    # Clean up intermediate files
    rm ${sam_file}
    """
}

process index {
    // index an aligned and sorted BAM
    input:
    path bam_file

    output:
    path "${bam_file}.bai"

    script:
    """
    samtools index ${bam_file}
    """
}
