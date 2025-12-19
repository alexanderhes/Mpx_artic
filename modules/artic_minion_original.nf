process ARTIC_MINION {
    tag "${sample_id}"
    publishDir "${params.output_dir}/${run_name}/minion/${barcode}_${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq), val(run_name), val(barcode)
    path bed_file
    path ref_file
    path artic_model_dir

    output:
    tuple val(run_name), path("${barcode}_${sample_id}.*"), emit: all_outputs
    tuple val(run_name), path("${barcode}_${sample_id}.consensus.fasta"), optional: true, emit: consensus
    tuple val(run_name), val(sample_id), val(barcode), path("${barcode}_${sample_id}.primertrimmed.rg.sorted.bam"), optional: true, emit: bam_file



    script:
    """

    artic minion \
        --model ${params.artic_model} \
        --model-dir ${artic_model_dir} \
        --bed ${bed_file} \
        --ref ${ref_file} \
        --read-file ${fastq} \
        ${barcode}_${sample_id}
    """
}