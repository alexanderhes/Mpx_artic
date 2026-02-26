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
include { USHER_DOWNLOAD }      from './modules/usher_download'
include { USHER_ALIGN }         from './modules/usher_align'
include { USHER_MSA }           from './modules/usher_msa'
include { USHER_MAKE_VCF }      from './modules/usher_make_vcf'
include { USHER_MASK_VCF }      from './modules/usher_mask_vcf'
include { USHER_PLACE }         from './modules/usher_place'
include { USHER_REPORT }        from './modules/usher_report'
include { VERSION_ARTIC }       from './modules/versions'
include { VERSION_NEXTCLADE }   from './modules/versions'
include { VERSION_USHER }       from './modules/versions'
include { VERSION_GOFASTA }     from './modules/versions'
include { VERSION_BEDTOOLS }    from './modules/versions'
include { VERSION_R }           from './modules/versions'
include { VERSION_R_APE }       from './modules/versions'
include { COLLECT_VERSIONS }    from './modules/versions'


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

    // -------------------------------------------------------------------------
    // UShER phylogenetic placement
    // -------------------------------------------------------------------------

    // Download latest global mpox tree and metadata
    USHER_DOWNLOAD()

    // Align combined consensus sequences against the reference
    USHER_ALIGN(
        COMBINE_FASTA.out.combined_fasta,
        file(params.usher_ref)
    )

    // Build padded MSA from SAM
    USHER_MSA(
        USHER_ALIGN.out.sam_file,
        file(params.usher_ref)
    )

    // Generate VCF from MSA
    USHER_MAKE_VCF(
        USHER_MSA.out.msa_fasta
    )

    // Mask problem sites from VCF
    USHER_MASK_VCF(
        USHER_MAKE_VCF.out.raw_vcf,
        file(params.usher_mask)
    )

    // Place samples in global tree and generate reports
    usher_place_input = USHER_MASK_VCF.out.masked_vcf
        .combine(USHER_DOWNLOAD.out.pb_file)
        .combine(USHER_DOWNLOAD.out.metadata)

    USHER_PLACE(usher_place_input)

    USHER_PLACE.out.final_report.view { "UShER report: \$it" }

    // Closest-neighbor R report
    usher_report_input = USHER_PLACE.out.optimized_nwk
        .join(COMBINE_FASTA.out.combined_fasta)
        .combine(USHER_DOWNLOAD.out.metadata)

    USHER_REPORT(usher_report_input, file(params.input_dir))

    USHER_REPORT.out.neighbor_report.view { "Closest neighbor report: \$it" }

    // Run R script - combines sequencing QC, Nextclade and UShER neighbor results
    R_COMBINE_RESULTS(
        combined_stats_input,
        depth_combined_input,
        file(params.input_dir),          // samplesheet
        nextclade_input.map { it[1] },   // path to {run_name}.csv
        USHER_REPORT.out.neighbor_report.map { it[1] }  // path to closest neighbor TSV
    )

    R_COMBINE_RESULTS.out.final_results.view { "Final results: \$it" }    // -------------------------------------------------------------------------
    // Capture tool versions for retrospective auditability
    // -------------------------------------------------------------------------
    run_name_ch = input_samples.map { it[2] }.first()

    VERSION_ARTIC()
    VERSION_NEXTCLADE()
    VERSION_USHER()
    VERSION_GOFASTA()
    VERSION_BEDTOOLS()
    VERSION_R()
    VERSION_R_APE()

    all_versions = VERSION_ARTIC.out
        .mix(
            VERSION_NEXTCLADE.out,
            VERSION_USHER.out,
            VERSION_GOFASTA.out,
            VERSION_BEDTOOLS.out,
            VERSION_R.out,
            VERSION_R_APE.out
        )
        .collect()

    COLLECT_VERSIONS(run_name_ch, all_versions)
}
