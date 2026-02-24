process USHER_ALIGN {
    container 'quay.io/biocontainers/artic:1.8.5--pyhdfd78af_0'

    input:
    tuple val(run_name), path(query_fasta)
    path ref_fasta

    output:
    tuple val(run_name), path('aligned.sam'), emit: sam_file

    script:
    """
    minimap2 -a -x asm20 --for-only --junc-bonus=0 --score-N=0 --secondary=no \\
        "${ref_fasta}" "${query_fasta}" > aligned.sam
    """
}
