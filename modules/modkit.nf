process modkit_pileup {
    // execute the modkit pileup command
    input:
    path bam_file
    path output_dir

    output:
    path "${output_dir}/modkit_pileup_output.bed" 

    script:
    """
    modkit pileup -t $task.cpus  ${bam_file} ${output_dir}/modkit_pileup_output.bed --only-tabs --filter-threshold ${params.filter_threshold_modkit}

    """
}
