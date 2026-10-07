process CONSENSUS_SUMMARY {
    publishDir "${params.output_dir}/${run_name}/", mode: 'copy'

    input:
    tuple val(run_name), path(nextclade_csv), path(neighbor_report)
    path sample_sheet

    output:
    tuple val(run_name), path("${run_name}_final_results.csv"), emit: final_results

    script:
    """
    Rscript ${projectDir}/bin/consensus_summary.R \
        ${nextclade_csv} \
        ${sample_sheet} \
        ${neighbor_report} \
        ${run_name}_final_results.csv
    """
}
