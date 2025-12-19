process GUPPYPLEX {
    tag "$sample_id"
    publishDir "${params.output_dir}/${run_name}/filtering", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_dir), val(run_name), val(barcode)

    output:
    tuple val(sample_id), path("${barcode}_${sample_id}_filtered.fastq"), val(run_name), val(barcode), emit: filtered_fastq

    script:
    """
    
    artic guppyplex \
        --directory "${fastq_dir}" \
        --min-length ${params.min_length} \
        --prefix "${barcode}_${sample_id}" \
        --output "${barcode}_${sample_id}_filtered.fastq"
    """
}