process COMBINE_FASTA {
    publishDir "${params.output_dir}/${run_name}/combined", mode: 'copy'

    input:
    tuple val(run_name), path(fasta_files), val(output_name)

    output:
    tuple val(run_name), path("${output_name}"), emit: combined_fasta

    script:
    """
    cat ${fasta_files.join(' ')} > ${output_name}
    """
}