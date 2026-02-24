process USHER_MASK_VCF {
    container 'quay.io/biocontainers/bedtools:2.31.1--h13024bc_3'

    input:
    tuple val(run_name), path(raw_vcf)
    path mask_bed

    output:
    tuple val(run_name), path('query.masked.vcf'), emit: masked_vcf

    script:
    """
    bedtools intersect -a "${raw_vcf}" -b "${mask_bed}" -v -header > query.masked.vcf
    """
}
