process NEXTCLADE {
    publishDir "${params.output_dir}/${run_name}/nextclade_output", mode: 'copy'
    container 'nextstrain/nextclade:3.9.1'


    input:
    tuple val(run_name), path(combined_consensus)

    output:
    tuple val(run_name), path("nextclade_output"), emit: nextclade_results
    tuple val(run_name), path("nextclade_output/${run_name}.csv"), optional: true, emit: nextclade_csv

    script:
    """
    # Download the latest nextclade dataset (all clades) fresh each run
    nextclade dataset get --name MPXV --output-dir nextclade_dataset
    
    nextclade run ${combined_consensus} \
        -D nextclade_dataset \
        -n ${run_name} \
        -O nextclade_output/

    echo "Nextclade analysis completed for ${run_name}"
    ls -l nextclade_output/
    """
}