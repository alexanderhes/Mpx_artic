// ─────────────────────────────────────────────────────────────────────────────
// Capture tool versions at runtime for retrospective auditability.
// One lightweight process per container; each emits CSV rows (tool,version,container).
// COLLECT_VERSIONS merges them into {run_name}_versions.csv in the output dir.
// ─────────────────────────────────────────────────────────────────────────────

process VERSION_ARTIC {
    container 'quay.io/biocontainers/artic:1.8.5--pyhdfd78af_0'

    output:
    path "versions_artic.csv"

    script:
    """
    VER=\$(artic --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+' | head -1 || echo "unknown")
    echo "artic,\${VER},quay.io/biocontainers/artic:1.8.5--pyhdfd78af_0" > versions_artic.csv
    """
}

process VERSION_NEXTCLADE {
    container 'nextstrain/nextclade:3.9.1'

    output:
    path "versions_nextclade.csv"

    script:
    """
    VER=\$(nextclade --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+' | head -1 || echo "unknown")
    echo "nextclade,\${VER},nextstrain/nextclade:3.9.1" > versions_nextclade.csv
    """
}

process VERSION_USHER {
    container 'quay.io/biocontainers/usher:0.6.6--hdd55de9_4'

    output:
    path "versions_usher.csv"

    script:
    """
    VER=\$(usher --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+' | head -1 || echo "unknown")
    echo "usher,\${VER},quay.io/biocontainers/usher:0.6.6--hdd55de9_4" > versions_usher.csv
    """
}

process VERSION_GOFASTA {
    container 'quay.io/biocontainers/gofasta:1.2.3--h9ee0642_0'

    output:
    path "versions_gofasta.csv"

    script:
    """
    VER=\$(gofasta --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+' | head -1 || echo "unknown")
    echo "gofasta,\${VER},quay.io/biocontainers/gofasta:1.2.3--h9ee0642_0" > versions_gofasta.csv
    """
}

process VERSION_BEDTOOLS {
    container 'quay.io/biocontainers/bedtools:2.31.1--h13024bc_3'

    output:
    path "versions_bedtools.csv"

    script:
    """
    VER=\$(bedtools --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+' | head -1 || echo "unknown")
    echo "bedtools,\${VER},quay.io/biocontainers/bedtools:2.31.1--h13024bc_3" > versions_bedtools.csv
    """
}

process VERSION_R {
    container 'rocker/tidyverse:4.4'

    output:
    path "versions_r.csv"

    script:
    """
    R_VER=\$(R --version 2>&1 | head -1 | grep -oP '\\d+\\.\\d+\\.\\d+' | head -1 || echo "unknown")
    echo "R,\${R_VER},rocker/tidyverse:4.4" > versions_r.csv
    Rscript -e '
        pkgs <- c("tidyverse","dplyr","readr","stringr","lubridate")
        for (p in pkgs) {
            cat(p, as.character(packageVersion(p)), "rocker/tidyverse:4.4", sep=",")
            cat("\\n")
        }
    ' 2>&1 >> versions_r.csv || true
    """
}

process VERSION_R_APE {
    container 'community.wave.seqera.io/library/r-ape_r-base_r-castor_r-lubridate_r-tidyverse:b29271c460ff37dc'

    output:
    path "versions_r_ape.csv"

    script:
    """
    R_VER=\$(R --version 2>&1 | head -1 | grep -oP '\\d+\\.\\d+\\.\\d+' | head -1 || echo "unknown")
    echo "R,\${R_VER},community.wave.seqera.io/library/r-ape_r-base_r-castor_r-lubridate_r-tidyverse:b29271c460ff37dc" > versions_r_ape.csv
    Rscript -e '
        pkgs <- c("ape","castor","lubridate","tidyverse")
        for (p in pkgs) {
            cat(p, as.character(packageVersion(p)), "community.wave.seqera.io/library/r-ape_r-base_r-castor_r-lubridate_r-tidyverse:b29271c460ff37dc", sep=",")
            cat("\\n")
        }
    ' 2>&1 >> versions_r_ape.csv || true
    """
}

process COLLECT_VERSIONS {
    publishDir "${params.output_dir}/${run_name}", mode: 'copy'

    input:
    val  run_name
    path version_files

    output:
    path "${run_name}_versions.csv"

    script:
    """
    {
        echo "tool,version,container,run_name,date"
        awk -v run="${run_name}" -v date="\$(date +'%Y-%m-%d')" '{print \$0","run","date}' ${version_files}
    } > ${run_name}_versions.csv
    """
}
