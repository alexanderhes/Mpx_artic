process R_COMBINE_RESULTS {
    publishDir "${params.output_dir}/${run_name}/", mode: 'copy'

    input:
    tuple val(run_name), path(combined_stats)
    tuple val(run_name), path(combined_depth_stats)
    path sample_sheet
    path nextclade_csv
    
    output:
    tuple val(run_name), path("${run_name}_final_results.csv"), emit: final_results


    script:
    """
    Rscript ${projectDir}/bin/R_summary.R \
        ${combined_stats} \
        ${sample_sheet} \
        ${nextclade_csv} \
        ${run_name}_final_results.csv \
        ${combined_depth_stats}
    """
}