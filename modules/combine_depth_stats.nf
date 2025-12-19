process COMBINE_DEPTH_STATS {
    debug true
    publishDir "${params.output_dir}/${run_name}/minion/depth_stats", mode: 'copy'

    input:
    // grouped by run_name like COMBINE_STATS does
    tuple val(run_name), val(barcode_sample_ids), path(depth_files)

    output:
    tuple val(run_name), path("${run_name}_combined_depth_stats.tsv"), emit: combined_depth_stats

    script:
    """
    echo "Debug: Combining depth stats for run ${run_name}"

    # Get header from first file
    first_file=\$(echo ${depth_files} | awk '{print \$1}')
    head -n 1 "\$first_file" > ${run_name}_combined_depth_stats.tsv

    # Append all rows (excluding headers)
    for file in ${depth_files}; do
        tail -n +2 "\$file" >> ${run_name}_combined_depth_stats.tsv
    done

    echo "Debug: Depth stats combined for run ${run_name}"
    """
}