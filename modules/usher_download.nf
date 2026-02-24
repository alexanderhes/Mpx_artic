process USHER_DOWNLOAD {
    container 'quay.io/biocontainers/artic:1.8.5--pyhdfd78af_0'

    output:
    path 'mpxv.latest.masked.pb',    emit: pb_file
    path 'mpxv.latest.metadata.tsv', emit: metadata

    script:
    def base_url = "https://hgdownload.gi.ucsc.edu/hubs/GCF/014/621/545/GCF_014621545.1/UShER_hMPXV"
    """
    curl -fsSL "${base_url}/mpxv.latest.masked.pb.gz"    | gunzip -c > mpxv.latest.masked.pb
    curl -fsSL "${base_url}/mpxv.latest.metadata.tsv.gz" | gunzip -c > mpxv.latest.metadata.tsv
    """
}
