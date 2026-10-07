#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import workflows
include { ARTIC }               from './workflows/artic'
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
    ARTIC()
    run_name_ch = ARTIC.out.run_name

    // -------------------------------------------------------------------------
    // Capture tool versions for retrospective auditability
    // -------------------------------------------------------------------------
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
