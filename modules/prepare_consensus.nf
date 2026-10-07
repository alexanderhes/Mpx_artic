process PREPARE_CONSENSUS {
    publishDir "${params.output_dir}/${run_name}/combined", mode: 'copy'

    input:
    tuple val(run_name), path(fasta_files, stageAs: 'input_fasta/*')
    path samplesheet
    path ref_fasta

    output:
    tuple val(run_name), path("${run_name}_combined_consensus.fasta"), emit: combined_fasta
    path "${run_name}_samplesheet_consensus.csv",                       emit: samplesheet
    path "${run_name}_input_check.tsv",                                 emit: input_check

    script:
    """
    python3 ${projectDir}/bin/prepare_consensus.py \\
        --samplesheet "${samplesheet}" \\
        --ref "${ref_fasta}" \\
        --run-name "${run_name}" \\
        --fasta input_fasta/*
    """
}
