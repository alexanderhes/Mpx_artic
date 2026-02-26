process USHER_MAKE_VCF {
    container 'quay.io/biocontainers/usher:0.6.6--hdd55de9_4'

    input:
    tuple val(run_name), path(msa_fasta)

    output:
    tuple val(run_name), path('raw.vcf'), emit: raw_vcf

    script:
    """
    faToVcf "${msa_fasta}" raw.vcf
    """
}
