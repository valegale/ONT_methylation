#!/usr/bin/env nextflow

process hello {
    input:
    val name

    output:
    stdout

    script:
    """
    echo 'Hello, ${name}!'
    """
}

workflow {
    hello(params.name ?: 'world')
}

