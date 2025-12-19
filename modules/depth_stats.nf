process DEPTH_STATS {
    debug true
    publishDir "${params.output_dir}/${run_name}/minion/depth_stats", mode: 'copy'

    input:
    // same tuple structure as READ_STATS, but BAM is the .sorted.bam
    tuple val(sample_id), path(raw_fastq), val(run_name), val(barcode), path(filtered_fastq), path(raw_bam)

    output:
    // per-sample TSV with unnormalized stats + average depth
    tuple val(sample_id), val(run_name), val(barcode), path("${barcode}_${sample_id}_depth_stats.tsv"), emit: depth_stats

    script:
    """
    echo "Debug: Starting DEPTH_STATS for sample ${sample_id}"

    # Count raw reads across all fastq.gz files in the directory
    raw_count=0
    for file in ${raw_fastq}/*.fastq.gz; do
        if [ -f "\$file" ]; then
            file_count=\$(zcat "\$file" | wc -l)
            file_count=\$((file_count / 4))
            raw_count=\$((raw_count + file_count))
        fi
    done

    # Count filtered reads
    filtered_count=\$(( \$(wc -l < ${filtered_fastq}) / 4 ))

    # Count mapped reads from UNNORMALIZED BAM
    mapped_count=\$(samtools view -F 2308 -c ${raw_bam})

    # Mapping percentage based on filtered reads
    mapped_percentage=\$(awk "BEGIN {print (\$mapped_count/\$filtered_count) * 100}")

    # Average depth across the reference using samtools depth
    avg_depth=\$(samtools depth -a ${raw_bam} | awk '{sum+=\$3; n++} END { if (n>0) print sum/n; else print 0 }')

    # Create output table
    echo -e "Sample\tRaw_Reads\tFiltered_Reads\tMapped_Reads_Unnormalized\tMapped_Percentage_Unnormalized\tAverageDepth" > ${barcode}_${sample_id}_depth_stats.tsv
    echo -e "${barcode}_${sample_id}\t\${raw_count}\t\${filtered_count}\t\${mapped_count}\t\${mapped_percentage}\t\${avg_depth}" >> ${barcode}_${sample_id}_depth_stats.tsv

    echo "Debug: DEPTH_STATS completed for sample ${sample_id}"
    """
}
