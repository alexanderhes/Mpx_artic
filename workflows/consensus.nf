// -------------------------------------------------------------------------
// Consensus mode: finished genomes (one FASTA per sample, file name = PrøveID)
// -> Nextclade typing and UShER placement in the global tree
// -------------------------------------------------------------------------

include { PREPARE_CONSENSUS }   from '../modules/prepare_consensus'
include { NEXTCLADE }           from '../modules/nextclade'
include { CONSENSUS_SUMMARY }   from '../modules/consensus_summary'
include { USHER_PLACEMENT }     from '../subworkflows/usher_placement'

workflow CONSENSUS {
    main:
    def samplesheet = file(params.input_dir, checkIfExists: true)

    // Run name comes from the samplesheet and must be a single value
    def run_names = samplesheet
        .splitCsv(header: true, sep: ';')
        .collect { row -> row.find { k, v -> k.replace('﻿', '').trim() == 'RunName' }?.value?.trim() }
        .unique()
    if (run_names.size() != 1 || !run_names[0]) {
        error "Consensus mode: samplesheet must have a RunName column with exactly one value, found: ${run_names}"
    }
    def run = run_names[0]

    fasta_files = channel
        .fromPath(params.fasta, checkIfExists: true)
        .collect()
        .map { files -> tuple(run, files) }

    // Validate inputs, rename headers to sample IDs and orient to the reference
    PREPARE_CONSENSUS(fasta_files, samplesheet, file(params.usher_ref))

    prepared_samplesheet = PREPARE_CONSENSUS.out.samplesheet

    // Clade/lineage assignment
    NEXTCLADE(PREPARE_CONSENSUS.out.combined_fasta)

    // UShER phylogenetic placement and closest-neighbor report
    USHER_PLACEMENT(PREPARE_CONSENSUS.out.combined_fasta, prepared_samplesheet)

    // Combine Nextclade and UShER neighbor results with the samplesheet
    CONSENSUS_SUMMARY(
        NEXTCLADE.out.nextclade_csv.join(USHER_PLACEMENT.out.neighbor_report),
        prepared_samplesheet
    )

    CONSENSUS_SUMMARY.out.final_results.view { "Final results: ${it[1]}" }

    emit:
    run_name = channel.value(run)
}
