process compute_statistics {
    // run an inhouse script that computes which bases are methylated: modified bases / total bases 
    publishDir  params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(bed_file), path(reference)

    output:
    path("${sample_id}/methylation_statistics")

    script:
    """
    python ${baseDir}/scripts/compute_statistics.py ${bed_file} ${reference} ${sample_id}/methylation_statistics
    """
}