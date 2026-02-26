process USHER_REPORT {
    container 'community.wave.seqera.io/library/r-ape_r-base_r-castor_r-lubridate_r-tidyverse:b29271c460ff37dc'
    publishDir "${params.output_dir}/${run_name}/usher", mode: 'copy'

    input:
    tuple val(run_name), path(optimized_nwk), path(combined_fasta), path(global_metadata)
    path samplesheet

    output:
    tuple val(run_name), path("${run_name}_closest_neighbor_report.tsv"), emit: neighbor_report

    script:
    """
    Rscript ${projectDir}/bin/usher_report.R \
        "${combined_fasta}" \
        "${samplesheet}" \
        "${optimized_nwk}" \
        "${global_metadata}" \
        "${run_name}_closest_neighbor_report.tsv"
    """
}
