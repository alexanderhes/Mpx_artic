## Quick start guide

A wrapper that runs a MPXV ARTIC Nextflow pipeline for whole‑genome analysis of Monkeypox virus (MPXV) data generated on Oxford Nanopore devices.

The pipeline is designed for the tiling amplicon scheme from Welkers et al. (https://www.protocols.io/view/monkeypox-virus-whole-genome-sequencing-using-comb-n2bvj6155lk5/v1). It uses the ARTIC tools (https://github.com/artic-network/fieldbioinformatics) to generate consensus genomes and Nextclade (https://github.com/nextstrain/nextclade) for clade/lineage assignment.

Additionally, the pipeline includes **phylogenetic placement** using UShER (Ultrafast Sample placement on Existing tRees), which places each sample into the global MPXV phylogenetic tree. This provides phylogenetic context including identification of the closest global neighbors, parsimony distance, and date ranges of related sequences.

### Pipeline features
- **Consensus generation**: ARTIC minion with amplicon normalisation and primer trimming
- **Quality metrics**: Read statistics, mapping depth, and genome coverage
- **Clade/Lineage assignment**: Nextclade analysis for MPXV classification
- **Phylogenetic placement**: UShER global tree integration with closest-neighbor reporting

### Key output files
- `{RunName}_final_results.csv` — Master summary table (QC, clade, phylogenetic distance)
- `{RunName}_full_tree_optimized.nwk` — Optimised phylogenetic tree (Newick format)
- `{RunName}_context_tree.nwk/.json` — Context subtree with closest global neighbors
- `{RunName}_closest_neighbor_report.tsv` — Per-sample phylogenetic neighbor summary

For detailed pipeline documentation, see [MANUAL.md](MANUAL.md).

#### How to run example 
```
screen -S mpxv_artic-TEST -d -m bash /home/ngs/ngs_scripts/mpxv_artic/wrapper.sh \
-r TEST \
-a MPX \
-y 2025
```
#### Options
* `-r, --run`
  Sequencing run name (e.g. `MPX012` or `TEST`). 

* `-a, --agens`
  Analysis group identifier (`MPX`).

* `-y, --year`
  Year of the sequencing run.

* `-m, --mode`
  `artic` (default) runs the full pipeline from raw FASTQ. `consensus` runs typing and phylogenetic placement on finished genomes (see below).

* `-f, --fasta-dir`
  N-drive folder with the consensus FASTA files. Required with `-m consensus`.

### Consensus mode (finished genomes)
Runs Nextclade typing and UShER placement in the global tree on genomes that are already assembled (e.g. SPAdes assemblies or genomes from other pipelines), skipping ARTIC.

- One FASTA file per sample, with **file name = PrøveID** (e.g. `SAMPLE001.fasta`; `.fasta`, `.fa`, `.fna`, optionally `.gz`). Headers are ignored and replaced by the PrøveID.
- Each file must contain exactly one sequence. Multi-contig assemblies must be scaffolded/merged first.
- Sequences that are reverse-complemented relative to the reference are flipped automatically.
- The run stops before analysis if a FASTA has no samplesheet row (or vice versa), an ID is duplicated, or SampleDate is not DD.MM.YYYY.

With the wrapper:
```
screen -S mpxv_consensus-RUN -d -m bash /home/ngs/ngs_scripts/mpxv_artic/wrapper.sh \
-r RUN \
-a MPX \
-y 2026 \
-m consensus \
-f /Virologi/NGS/<folder>/<RUN>
```
`-f` is the N-drive folder (on the same share as the other wrapper paths) that holds the FASTA files; the files are copied locally before the run.
Directly with Nextflow:
```
nextflow run alexanderhes/Mpx_artic -r main \
    --fasta '/path/to/genomes/*.{fa,fasta,fna}' \
    --input_dir /path/to/RUN_samplesheet.csv
```

### Troubleshooting
Check status file 
```
cat ~/mpx_${RUN}_status.txt
```
View logs
```
tail -f ~/mpx_wrapper.log
tail -f ~/mpx_wrapper_error.log
```

## Samplesheet information

The pipeline requires a `{RunName}_samplesheet.csv` as input which must be ";" separated (default separator when saving CSV files from Norwegian Excel versions).

### Required columns
- **PrøveID** — Sample identification (e.g., `sample1`, `sample2`). In consensus mode this must match the FASTA file name.
- **barcode** — Sequencing barcode (e.g., `Barcode65`, `Barcode66`). Not required in consensus mode.
- **RunName** — Name of the sequencing run (must match the pipeline run name, e.g., `TEST`, `MPX012`)
- **SampleDate** — Collection date in DD.MM.YYYY format (e.g., `10.07.2022`) — **Required for phylogenetic analysis. Script will fail if missing.**

### Optional columns
- **PrøveNr** — Internal sample number
- **SequenceID** — Alternative sequence identifier
- **Ct** — PCR Ct value (if available)
- **Fortynning** — Sample dilution information
- **Kommentar** — Comments or notes
- **RunDate** — Date when sequencing was performed

### Samplesheet example 
```
PrøveID;SequenceID;Ct;Fortynning;SampleDate;Kommentar;RunDate;RunName;barcode
sample1;SEQ001;19;;10.07.2022;kommentar1;12.02.2025;TEST;Barcode65
sample3;SEQ002;20;;10.07.2022;kommentar3;12.02.2025;TEST;Barcode66
sample5;SEQ003;;;15.07.2022;;13.02.2025;TEST;Barcode67
```

**Note:** SampleDate is used to contextualize results alongside closest global neighbors in the UShER phylogenetic analysis. 
