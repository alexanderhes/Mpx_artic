// ─────────────────────────────────────────────────────────────────────────────
// Capture tool versions at runtime for retrospective auditability.
// One lightweight process per container; each emits a small text snippet.
// COLLECT_VERSIONS merges them into {run_name}_versions.txt in the output dir.
// ─────────────────────────────────────────────────────────────────────────────

process VERSION_ARTIC {
    container 'quay.io/biocontainers/artic:1.8.5--pyhdfd78af_0'

    output:
    path "versions_artic.txt"

    script:
    """
    {
        echo "=== ARTIC fieldbioinformatics ==="
        artic --version 2>&1 || true
        echo ""
        echo "=== minimap2 ==="
        minimap2 --version 2>&1 || true
        echo ""
        echo "=== samtools ==="
        samtools --version 2>&1 | head -1 || true
        echo ""
        echo "=== medaka ==="
        medaka --version 2>&1 || true
    } > versions_artic.txt
    """
}

process VERSION_NEXTCLADE {
    container 'nextstrain/nextclade:3.9.1'

    output:
    path "versions_nextclade.txt"

    script:
    """
    {
        echo "=== Nextclade ==="
        nextclade --version 2>&1 || true
    } > versions_nextclade.txt
    """
}

process VERSION_USHER {
    container 'quay.io/biocontainers/usher:0.6.6--hdd55de9_4'

    output:
    path "versions_usher.txt"

    script:
    """
    {
        echo "=== UShER ==="
        usher --version 2>&1 || true
        echo ""
        echo "=== matOptimize ==="
        matOptimize --version 2>&1 || true
        echo ""
        echo "=== matUtils ==="
        matUtils --version 2>&1 || true
        echo ""
        echo "=== faToVcf ==="
        faToVcf 2>&1 | head -2 || true
    } > versions_usher.txt
    """
}

process VERSION_GOFASTA {
    container 'quay.io/biocontainers/gofasta:1.2.3--h9ee0642_0'

    output:
    path "versions_gofasta.txt"

    script:
    """
    {
        echo "=== gofasta ==="
        gofasta --version 2>&1 || true
    } > versions_gofasta.txt
    """
}

process VERSION_BEDTOOLS {
    container 'quay.io/biocontainers/bedtools:2.31.1--h13024bc_3'

    output:
    path "versions_bedtools.txt"

    script:
    """
    {
        echo "=== bedtools ==="
        bedtools --version 2>&1 || true
    } > versions_bedtools.txt
    """
}

process VERSION_R {
    container 'rocker/tidyverse:4.4'

    output:
    path "versions_r.txt"

    script:
    """
    {
        echo "=== R ==="
        R --version 2>&1 | head -1 || true
        echo ""
        echo "=== R packages ==="
        Rscript -e 'pkgs <- c("tidyverse","dplyr","readr","stringr","lubridate"); for(p in pkgs) cat(p, as.character(packageVersion(p)), "\\n")' 2>&1 || true
    } > versions_r.txt
    """
}

process VERSION_R_APE {
    container 'community.wave.seqera.io/library/r-ape_r-base_r-castor_r-lubridate_r-tidyverse:b29271c460ff37dc'

    output:
    path "versions_r_ape.txt"

    script:
    """
    {
        echo "=== R (UShER report container) ==="
        R --version 2>&1 | head -1 || true
        echo ""
        echo "=== R packages ==="
        Rscript -e 'pkgs <- c("ape","castor","lubridate","tidyverse"); for(p in pkgs) cat(p, as.character(packageVersion(p)), "\\n")' 2>&1 || true
    } > versions_r_ape.txt
    """
}

process COLLECT_VERSIONS {
    publishDir "${params.output_dir}/${run_name}", mode: 'copy'

    input:
    val  run_name
    path version_files

    output:
    path "${run_name}_versions.txt"

    script:
    """
    {
        echo "MPXV ARTIC Pipeline — Tool Versions"
        echo "Run: ${run_name}"
        echo "Date: \$(date +'%Y-%m-%d %H:%M:%S')"
        echo "======================================="
        echo ""
        for f in ${version_files}; do
            cat "\$f"
            echo ""
        done
    } > ${run_name}_versions.txt
    """
}
