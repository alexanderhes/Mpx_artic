// -------------------------------------------------------------------------
// UShER phylogenetic placement
// Shared by the ARTIC and consensus workflows.
// -------------------------------------------------------------------------

include { USHER_DOWNLOAD }      from '../modules/usher_download'
include { USHER_ALIGN }         from '../modules/usher_align'
include { USHER_MSA }           from '../modules/usher_msa'
include { USHER_MAKE_VCF }      from '../modules/usher_make_vcf'
include { USHER_MASK_VCF }      from '../modules/usher_mask_vcf'
include { USHER_PLACE }         from '../modules/usher_place'
include { USHER_REPORT }        from '../modules/usher_report'

workflow USHER_PLACEMENT {
    take:
    combined_fasta  // (run_name, combined consensus fasta)
    samplesheet     // path to samplesheet (contains SampleDate)

    main:
    // Download latest global mpox tree and metadata
    USHER_DOWNLOAD()

    // Align combined consensus sequences against the reference
    USHER_ALIGN(
        combined_fasta,
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
        .join(combined_fasta)
        .combine(USHER_DOWNLOAD.out.metadata)

    USHER_REPORT(usher_report_input, samplesheet)

    USHER_REPORT.out.neighbor_report.view { "Closest neighbor report: \$it" }

    emit:
    neighbor_report = USHER_REPORT.out.neighbor_report  // (run_name, closest_neighbor_report.tsv)
    optimized_nwk   = USHER_PLACE.out.optimized_nwk
    final_report    = USHER_PLACE.out.final_report
}
