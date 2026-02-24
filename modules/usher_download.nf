process USHER_DOWNLOAD {
    container 'quay.io/biocontainers/usher:0.6.6--hdd55de9_4'

    output:
    path 'mpxv.latest.masked.pb',    emit: pb_file
    path 'mpxv.latest.metadata.tsv', emit: metadata

    script:
    def base_url = "https://hgdownload.gi.ucsc.edu/hubs/GCF/014/621/545/GCF_014621545.1/UShER_hMPXV"
    """
    python3 - << 'PYEOF'
import urllib.request, gzip, shutil

base = "${base_url}"
files = {
    "mpxv.latest.masked.pb.gz":    "mpxv.latest.masked.pb",
    "mpxv.latest.metadata.tsv.gz": "mpxv.latest.metadata.tsv",
}
for gz_name, out_name in files.items():
    url = base + "/" + gz_name
    print(f"Downloading {url}", flush=True)
    with urllib.request.urlopen(url) as response, \
         gzip.open(response) as gz_fh, \
         open(out_name, 'wb') as out_fh:
        shutil.copyfileobj(gz_fh, out_fh)
    print(f"Written: {out_name}", flush=True)
PYEOF
    """
}
