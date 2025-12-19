process COMBINE_STATS {
    publishDir "${params.output_dir}/${run_name}/minion/read_stats", mode: 'copy'

    input:
    tuple val(run_name), val(barcode_sample_ids), path(stats_files)

    output:
    tuple val(run_name), path("${run_name}_combined_read_stats.tsv"), emit: combined_stats

    script:
    """
    # Get header from first file
    head -n 1 ${stats_files[0]} > ${run_name}_combined_read_stats.tsv
    # Combine all files, excluding headers
    for file in ${stats_files}; do
        tail -n +2 \$file >> ${run_name}_combined_read_stats.tsv
    done
    """
}