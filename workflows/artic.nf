// -------------------------------------------------------------------------
// ARTIC mode: raw ONT FASTQ -> consensus -> QC, Nextclade and UShER placement
// -------------------------------------------------------------------------

include { GUPPYPLEX } from '../modules/guppyplex'
include { ARTIC_MINION } from '../modules/artic_minion'
include { COMBINE_FASTA } from '../modules/combine_fasta'
include { NEXTCLADE } from '../modules/nextclade'
include { READ_STATS } from '../modules/read_stats'
include { COMBINE_STATS } from '../modules/combine_stats'
include { R_COMBINE_RESULTS } from '../modules/combine_results'
include { DEPTH_STATS }         from '../modules/depth_stats'
include { COMBINE_DEPTH_STATS } from '../modules/combine_depth_stats'
include { USHER_PLACEMENT }     from '../subworkflows/usher_placement'

workflow ARTIC {
    main:
    // Read the sample sheet
    channel
        .fromPath(params.input_dir)
        .splitCsv(header:true, sep:';')
        .map { row ->
            def result = tuple(row.sample_id, file(row.fastq), row.RunName, row.barcode)
            if (params.debug) println "Debug: Mapped row: ${result}"
            return result
        }
        .set { input_samples }

    // Run Guppyplex
    GUPPYPLEX(input_samples)

    //Run Artic Minion
    ARTIC_MINION(
    GUPPYPLEX.out.filtered_fastq,
    file(params.bed_file),
    file(params.ref_file),
    file(params.artic_model_dir)
    )

    read_stats_input = input_samples
        .join(GUPPYPLEX.out.filtered_fastq, by: [0, 2, 3])
        .map { sample_id, run_name, barcode, raw_fastq, filtered_fastq ->
            [run_name, sample_id, barcode, raw_fastq, filtered_fastq]
        }
        .combine(ARTIC_MINION.out.bam_file, by: [0, 1, 2])
        .map { run_name, sample_id, barcode, raw_fastq, filtered_fastq, bam_file ->
            tuple(sample_id, raw_fastq, run_name, barcode, filtered_fastq, bam_file)
        }

    // Run READ_STATS
    READ_STATS(read_stats_input)

    //Create a new input channel for depth stats and depth from unormalized BAM files
    depth_stats_input = input_samples
    .join(GUPPYPLEX.out.filtered_fastq, by: [0, 2, 3])
    .map { sample_id, run_name, barcode, raw_fastq, filtered_fastq ->
        [run_name, sample_id, barcode, raw_fastq, filtered_fastq]
    }
    .combine(ARTIC_MINION.out.raw_bam, by: [0, 1, 2])
    .map { run_name, sample_id, barcode, raw_fastq, filtered_fastq, raw_bam ->
        tuple(sample_id, raw_fastq, run_name, barcode, filtered_fastq, raw_bam)
    }

    DEPTH_STATS(depth_stats_input)


    // Collect all consensus sequences

    consensus_sequences = ARTIC_MINION.out.consensus
        .map { run_name, fasta -> [run_name, fasta] }
        .groupTuple()

    combined_consensus = consensus_sequences.map { run_name, fastas ->
        [run_name, fastas.flatten(), "${run_name}_combined_consensus.fasta"]
    }

    // Run the new module with the combined consensus
    COMBINE_FASTA(combined_consensus)

    // Run Nextclade on combined consenses
    NEXTCLADE(COMBINE_FASTA.out.combined_fasta)


    // Collect all read stats
    all_read_stats = READ_STATS.out.read_stats
        .map { sample_id, run_name, barcode, stats_file ->
            [run_name, "${barcode}_${sample_id}", stats_file]
        }
        .groupTuple()

    // Combine all read stats into a single file
    COMBINE_STATS(all_read_stats)

    // Collect all depth and unormalized read stats
    all_depth_stats = DEPTH_STATS.out.depth_stats
        .map { sample_id, run_name, barcode, depth_stats ->
            [run_name, "${barcode}_${sample_id}", depth_stats]
        }
        .groupTuple()

    // Combine all depth and unormalized read stats into a single file
    COMBINE_DEPTH_STATS(all_depth_stats)

    // Prepare inputs for R_COMBINE_RESULTS
    combined_stats_input  = COMBINE_STATS.out.combined_stats              // (run_name, combined_read_stats.tsv)
    depth_combined_input  = COMBINE_DEPTH_STATS.out.combined_depth_stats  // (run_name, combined_depth_stats.tsv)
    nextclade_input       = NEXTCLADE.out.nextclade_csv                   // (run_name, nextclade_csv)

    // Optional sanity checks
    combined_stats_input.ifEmpty  { error "Combined stats input is empty" }
    depth_combined_input.ifEmpty  { error "Depth combined input is empty" }
    nextclade_input.ifEmpty       { error "Nextclade input is empty" }

    // UShER phylogenetic placement and closest-neighbor report
    USHER_PLACEMENT(
        COMBINE_FASTA.out.combined_fasta,
        file(params.input_dir)
    )

    // Run R script - combines sequencing QC, Nextclade and UShER neighbor results
    R_COMBINE_RESULTS(
        combined_stats_input,
        depth_combined_input,
        file(params.input_dir),          // samplesheet
        nextclade_input.map { it[1] },   // path to {run_name}.csv
        USHER_PLACEMENT.out.neighbor_report.map { it[1] }  // path to closest neighbor TSV
    )

    R_COMBINE_RESULTS.out.final_results.view { "Final results: \$it" }

    emit:
    run_name = input_samples.map { it[2] }.first()
}
