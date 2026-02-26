process USHER_PLACE {
    container 'quay.io/biocontainers/usher:0.6.6--hdd55de9_4'
    publishDir "${params.output_dir}/${run_name}/usher", mode: 'copy'

    cpus 8

    input:
    tuple val(run_name), path(masked_vcf), path(pb_file), path(metadata)

    output:
    tuple val(run_name), path("${run_name}_full_tree_global.nwk"),      emit: global_nwk
    tuple val(run_name), path("${run_name}_full_tree_global.json"),     emit: global_json
    tuple val(run_name), path("${run_name}_full_tree_optimized.nwk"),   emit: optimized_nwk
    tuple val(run_name), path("${run_name}_full_tree_optimized.json"),  emit: optimized_json
    tuple val(run_name), path("${run_name}_context_tree.nwk"),          emit: context_nwk
    tuple val(run_name), path("${run_name}_context_tree.json"),         emit: context_json
    tuple val(run_name), path("${run_name}_final_report.tsv"),          emit: final_report
    tuple val(run_name), path("${run_name}_final_report_10SNPs.tsv"),   emit: report_10snps

    script:
    """
    mkdir -p usher_results

    # Extract placed sample names from VCF column headers
    grep -m1 "^#CHROM" "${masked_vcf}" | cut -f10- | tr '\\t' '\\n' > sample_list.txt

    # -------------------------------------------------------------------------
    # Step 5: UShER global placement
    # -------------------------------------------------------------------------
    usher \\
        -i "${pb_file}" \\
        -v "${masked_vcf}" \\
        -k 0 \\
        -o usher_results/placed.pb \\
        -d usher_results/

    # Export global tree - Newick
    matUtils extract \\
        -i usher_results/placed.pb \\
        -t ${run_name}_full_tree_global.nwk

    # Export global tree - JSON with metadata
    matUtils extract \\
        -i usher_results/placed.pb \\
        -M "${metadata}" \\
        -j ${run_name}_full_tree_global.json

    # -------------------------------------------------------------------------
    # Step 5.5: Tree optimisation
    # -------------------------------------------------------------------------
    matOptimize \\
        -i usher_results/placed.pb \\
        -o usher_results/placed_optimized.pb \\
        -r 4 \\
        -T ${task.cpus}

    # Export optimized tree - Newick
    matUtils extract \\
        -i usher_results/placed_optimized.pb \\
        -t ${run_name}_full_tree_optimized.nwk

    # Export optimized tree - JSON with metadata
    matUtils extract \\
        -i usher_results/placed_optimized.pb \\
        -M "${metadata}" \\
        -j ${run_name}_full_tree_optimized.json

    # -------------------------------------------------------------------------
    # Step 6: Context reports
    # -------------------------------------------------------------------------

    # Closest relative report (absolute nearest)
    matUtils extract \\
        -i usher_results/placed_optimized.pb \\
        -s sample_list.txt \\
        -V ${run_name}_final_report.tsv

    # Closest relatives within 10 SNPs
    matUtils extract \\
        -i usher_results/placed_optimized.pb \\
        -s sample_list.txt \\
        --within-distance ${run_name}_final_report_10SNPs.tsv \\
        --distance-threshold 10

    # Context subtree - JSON coloured with metadata
    matUtils extract \\
        -i usher_results/placed_optimized.pb \\
        -s sample_list.txt \\
        -Y 20 \\
        -M "${metadata}" \\
        -j ${run_name}_context_tree.json

    # Context subtree - Newick (for R/ggtree)
    matUtils extract \\
        -i usher_results/placed_optimized.pb \\
        -s sample_list.txt \\
        -Y 20 \\
        -t ${run_name}_context_tree.nwk
    """
}
