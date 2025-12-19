#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { GUPPYPLEX } from './modules/guppyplex'
include { ARTIC_MINION } from './modules/artic_minion'
include { COMBINE_FASTA } from './modules/combine_fasta'
include { NEXTCLADE } from './modules/nextclade'
include { READ_STATS } from './modules/read_stats'
include { COMBINE_STATS } from './modules/combine_stats'
include { R_COMBINE_RESULTS } from './modules/combine_results'
include { DEPTH_STATS }         from './modules/depth_stats'
include { COMBINE_DEPTH_STATS } from './modules/combine_depth_stats'


// Define parameters
params.input_dir = "$projectDir/samplesheet_mpox_test.csv"

// Define the main workflow
workflow {
    // Read the sample sheet
    channel
        .fromPath(params.input_dir)
        .splitCsv(header:true, sep:';')
        .map { row -> 
            def result = tuple(row.sample_id, file(row.fastq), row.RunName, row.barcode)
            println "Debug: Mapped row: ${result}"
            return result
        }
        .set { input_samples }

    // Print debug information
    input_samples.view { sample_id, fastq_dir, run_name, barcode ->
        "Debug: sample_id: ${sample_id}, fastq_dir: ${fastq_dir.toString()}, run_name: ${run_name}, barcode: ${barcode}"
    }
    // Run Guppyplex
    GUPPYPLEX(input_samples)

    GUPPYPLEX.out.filtered_fastq.view{ "Debug: Guppy FastQ: $it" }

    //Run Artic Minion
    ARTIC_MINION(
    GUPPYPLEX.out.filtered_fastq,
    file(params.bed_file),
    file(params.ref_file),
    file(params.artic_model_dir)
    )
    
    input_samples.view { "Debug: input_samples: $it" }
    ARTIC_MINION.out.bam_file.view{ "Debug: BAM file: $it" }
    ARTIC_MINION.out.consensus.view { "Raw consensus: $it" }


    read_stats_input = input_samples
        .join(GUPPYPLEX.out.filtered_fastq, by: [0, 2, 3])
        .view { "Debug: After first join: $it" }
        .map { sample_id, run_name, barcode, raw_fastq, filtered_fastq ->
            [run_name, sample_id, barcode, raw_fastq, filtered_fastq]
        }
        .view { "Debug: Before combine: $it" }
        .combine(ARTIC_MINION.out.bam_file, by: [0, 1, 2])
        .view { "Debug: After combine: $it" }
        .map { run_name, sample_id, barcode, raw_fastq, filtered_fastq, bam_file ->
            tuple(sample_id, raw_fastq, run_name, barcode, filtered_fastq, bam_file)
        }
        .view { "Debug: read_stats_input: $it" }

    // Run READ_STATS
    READ_STATS(read_stats_input)
    READ_STATS.out.read_stats.view{ "Read stats file: $it"}

    //Create a new input channel for depth stats and depth from unormalized BAM files
    depth_stats_input = input_samples
    .join(GUPPYPLEX.out.filtered_fastq, by: [0, 2, 3])
    .view { "Debug: DEPTH After first join: $it" }
    .map { sample_id, run_name, barcode, raw_fastq, filtered_fastq ->
        [run_name, sample_id, barcode, raw_fastq, filtered_fastq]
    }
    .view { "Debug: DEPTH Before combine: $it" }
    .combine(ARTIC_MINION.out.raw_bam, by: [0, 1, 2])
    .view { "Debug: DEPTH After combine: $it" }
    .map { run_name, sample_id, barcode, raw_fastq, filtered_fastq, raw_bam ->
        tuple(sample_id, raw_fastq, run_name, barcode, filtered_fastq, raw_bam)
    }
    .view { "Debug: depth_stats_input: $it" }

    DEPTH_STATS(depth_stats_input)
    DEPTH_STATS.out.depth_stats.view { "Depth stats file: $it" }


    // Collect all consensus sequences

    consensus_sequences = ARTIC_MINION.out.consensus
        .map { run_name, fasta -> [run_name, fasta] }
        .groupTuple()
        .view { "Grouped consensus: $it" }

    combined_consensus = consensus_sequences.map { run_name, fastas -> 
        [run_name, fastas.flatten(), "${run_name}_combined_consensus.fasta"]
    }.view { "Combined consensus: $it" }

    // Run the new module with the combined consensus
    COMBINE_FASTA(combined_consensus)

    COMBINE_FASTA.out.combined_fasta.view()

    // Run Nextclade on combined consenses
    NEXTCLADE(COMBINE_FASTA.out.combined_fasta, file(params.nextclade_dataset, checkIfExists: true))

    NEXTCLADE.out.nextclade_csv.view { "Debug: NEXTCLADE output: $it" }


    // Collect all read stats
    all_read_stats = READ_STATS.out.read_stats
        .map { sample_id, run_name, barcode, stats_file -> 
            [run_name, "${barcode}_${sample_id}", stats_file]
        }
        .groupTuple()
        .view { "Debug: Grouped read stats: $it" }

    // Combine all read stats into a single file
    COMBINE_STATS(all_read_stats)
    COMBINE_STATS.out.combined_stats.view { "Combined stats: $it" }

    // Collect all depth and unormalized read stats
    all_depth_stats = DEPTH_STATS.out.depth_stats
        .map { sample_id, run_name, barcode, depth_stats -> 
            [run_name, "${barcode}_${sample_id}", depth_stats]
        }
        .groupTuple()
        .view { "Debug: Grouped read stats: $it" }

    // Combine all depth and unormalized read stats into a single file
    COMBINE_DEPTH_STATS(all_depth_stats)
    COMBINE_DEPTH_STATS.out.combined_depth_stats.view { "Depth stats: $it" }

    // Prepare inputs for R_COMBINE_RESULTS
    combined_stats_input  = COMBINE_STATS.out.combined_stats              // (run_name, combined_read_stats.tsv)
    depth_combined_input  = COMBINE_DEPTH_STATS.out.combined_depth_stats  // (run_name, combined_depth_stats.tsv)
    nextclade_input       = NEXTCLADE.out.nextclade_csv                   // (run_name, nextclade_csv)

    combined_stats_input.view  { "Debug: COMBINED_STATS input: $it" }
    depth_combined_input.view  { "Debug: DEPTH_COMBINED input: $it" }
    nextclade_input.view       { "Debug: NEXTCLADE input: $it" }

    // Optional sanity checks
    combined_stats_input.ifEmpty  { error "Combined stats input is empty" }
    depth_combined_input.ifEmpty  { error "Depth combined input is empty" }
    nextclade_input.ifEmpty       { error "Nextclade input is empty" }

    // Run R script to combine results
    R_COMBINE_RESULTS(
        combined_stats_input,
        depth_combined_input,
        file(params.input_dir),          // samplesheet
        nextclade_input.map { it[1] }    // path to {run_name}.csv
    )

    // View the final results
    R_COMBINE_RESULTS.out.final_results.view { "Final results: $it" }

}
