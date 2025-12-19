process READ_STATS {
    debug true
    publishDir "${params.output_dir}/${run_name}/minion/read_stats", mode: 'copy'
    
    input:
    tuple val(sample_id), path(raw_fastq), val(run_name), val(barcode), path(filtered_fastq), path(bam_file)

    output:
    tuple val(sample_id), val(run_name), val(barcode), path("${barcode}_${sample_id}_read_stats.tsv"), emit: read_stats

    script:
    """
    echo "Debug: Starting READ_STATS for sample ${sample_id}"
    # Count raw reads across all fastq.gz files in the directory
    raw_count=0
    for file in ${raw_fastq}/*.fastq.gz; do
        file_count=\$(zcat \$file | echo \$((`wc -l`/4)))
        raw_count=\$((raw_count + file_count))
    done

    # Count filtered reads
    filtered_count=\$(echo \$((`wc -l < ${filtered_fastq}`/4)))

    # Count mapped reads and calculate percentage
    mapped_count=\$(samtools view -F 2308 -c ${bam_file})
    mapped_percentage=\$(awk "BEGIN {print (\$mapped_count/\$filtered_count) * 100}")

    # Create output table
    echo -e "Sample\tRaw_Reads\tFiltered_Reads\tMapped_Reads\tMapped_Percentage" > ${barcode}_${sample_id}_read_stats.tsv
    echo -e "${barcode}_${sample_id}\t\${raw_count}\t\${filtered_count}\t\${mapped_count}\t\${mapped_percentage}" >> ${barcode}_${sample_id}_read_stats.tsv
    """
}