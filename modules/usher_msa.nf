process USHER_MSA {
    container 'quay.io/biocontainers/gofasta:1.2.3--h9ee0642_0'

    input:
    tuple val(run_name), path(sam_file)
    path ref_fasta

    output:
    tuple val(run_name), path('aligned.msa.fasta'), emit: msa_fasta

    script:
    """
    gofasta sam toMultiAlign \\
        -s "${sam_file}" \\
        --reference "${ref_fasta}" \\
        --pad \\
        -o temp_queries.fasta

    cat "${ref_fasta}" temp_queries.fasta > aligned.msa.fasta
    rm temp_queries.fasta
    """
}
