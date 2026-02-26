# MPXV ARTIC Pipeline — Manual

## Overview

This pipeline performs whole-genome analysis of Monkeypox virus (MPXV) sequencing data generated on Oxford Nanopore devices using the Welkers et al. tiling amplicon scheme. It produces consensus genomes, quality metrics, clade/lineage assignments and phylogenetic placement into the global MPXV tree.

---

## Pipeline steps

### Step 1 — Read filtering (`GUPPYPLEX`)
Raw ONT reads from each barcode directory are length-filtered using `artic guppyplex`. Reads shorter than `params.min_length` (default: 100 bp) are discarded. 

**Input:** Per-sample raw FASTQ directories (from samplesheet)  
**Output:** One filtered FASTQ file per sample (`{barcode}_{sample_id}_filtered.fastq`)

---

### Step 2 — Consensus genome assembly (`ARTIC_MINION`)
Filtered reads are assembled into consensus genomes using `artic minion`. The tool performs:
- Read mapping to the reference genome (`NC_063383_masked.fasta`)
- Amplicon normalisation and primer trimming
- Variant calling using the specified basecalling model
- Consensus FASTA generation with low-coverage positions masked as N

**Input:** Filtered FASTQ, primer BED file, reference FASTA, basecalling model  
**Output:** Per-sample consensus FASTA, primer-trimmed BAM, raw (un-normalised) BAM, and all intermediate artic minion files

---

### Step 3 — Read and mapping statistics (`READ_STATS`, `DEPTH_STATS`)
Two parallel statistics modules run on each sample:

- **READ_STATS**: Counts raw reads, filtered reads, and reads mapped to the reference from the primer-trimmed BAM
- **DEPTH_STATS**: Calculates average depth and mapping statistics from the un-normalised BAM, reflecting true sequencing depth before amplicon normalisation

**Input:** Raw FASTQ, filtered FASTQ, BAM files  
**Output:** Per-sample TSV statistics files

---

### Step 4 — Combine consensus sequences (`COMBINE_FASTA`)
All per-sample consensus FASTAs from a run are concatenated into a single multi-FASTA file. This is used as input for both Nextclade and the UShER phylogenetic pipeline.

**Input:** All per-sample consensus FASTAs for the run  
**Output:** `{RunName}_combined_consensus.fasta`

---

### Step 5 — Clade and lineage assignment (`NEXTCLADE`)
The combined consensus FASTA is analysed with Nextclade using the local mpox dataset. Nextclade assigns each sample to an MPXV clade and lineage, and calculates genome coverage.

**Input:** Combined consensus FASTA, Nextclade dataset  
**Output:** Nextclade results folder including a CSV with per-sample assignments

---

### Step 6 — Download global MPXV tree (`USHER_DOWNLOAD`)
The latest global MPXV masked phylogenetic tree and metadata are downloaded fresh from the UCSC UShER hub at the start of each run, ensuring samples are placed in the most up-to-date global context.

**Source:** `https://hgdownload.gi.ucsc.edu/hubs/GCF/014/621/545/GCF_014621545.1/UShER_hMPXV/`  
**Output:** `mpxv.latest.masked.pb` (protobuf tree), `mpxv.latest.metadata.tsv`

---

### Step 7 — Alignment to reference (`USHER_ALIGN`)
The combined consensus FASTA is aligned to the static MPXV reference genome (`assets/UShER/NC_063383.fasta`) using `minimap2` with assembly-to-reference parameters

**Input:** Combined consensus FASTA, reference FASTA  
**Output:** SAM alignment file

---

### Step 8 — Multiple sequence alignment (`USHER_MSA`)
The SAM alignment is converted to a padded multiple sequence alignment (MSA) FASTA using `gofasta sam toMultiAlign`. The reference sequence is prepended, which is required for VCF generation in the next step.

**Input:** SAM file, reference FASTA  
**Output:** `aligned.msa.fasta`

---

### Step 9 — VCF generation (`USHER_MAKE_VCF`)
A VCF file of variant positions is generated from the MSA using `faToVcf`.

**Input:** MSA FASTA  
**Output:** `raw.vcf`

---

### Step 10 — VCF masking (`USHER_MASK_VCF`)
Known problematic sites (recombination hotspots, homopolymers, etc.) are removed from the VCF using `bedtools intersect` with the static mask file (`assets/UShER/mask_cladeII.bed`). This prevents artefactual placements in the phylogeny.

**Input:** Raw VCF, mask BED  
**Output:** `query.masked.vcf`

---

### Step 11 — Phylogenetic placement and tree export (`USHER_PLACE`)
Samples are placed into the global MPXV tree using UShER (Ultrafast Sample placement on Existing tRees):

1. **Global placement** (`usher`): Samples are placed using a global algorithm minimising parsimony score
2. **Tree optimisation** (`matOptimize -r 4`): The tree topology is further optimised with 4 rounds of stochastic nearest-neighbour interchange
3. **Tree export** (`matUtils extract`): Both global and optimised trees are exported in Newick and JSON formats
4. **Context reports**: Per-sample closest-relative reports are generated — one listing the 20 nearest global neighbours per local sample, one listing all neighbours within 10 SNPs

**Input:** Masked VCF, global `.pb` tree, global metadata TSV  
**Output:** See output files section below

---

### Step 12 — Closest-neighbour summary report (`USHER_REPORT`)
An R script calculates the phylogenetic distance from each local sample to all global samples in the optimised tree and produces a human-readable summary. For each local sample it reports:
- Number of equally close global neighbours (ties)
- Parsimony distance to the closest neighbour(s)
- Names of all closest neighbours
- Countries represented among closest neighbours
- Date range of closest neighbours

**Input:** Combined consensus FASTA (for tip label extraction), samplesheet (for sample dates), optimised Newick tree, global metadata TSV  
**Output:** `{RunName}_closest_neighbor_report.tsv`

---

### Step 13 — Final summary report (`R_COMBINE_RESULTS`)
An R script merges all per-run data into a single result table. This is the final step before version collection, combining all molecular and phylogenetic data.

Data merged:
- Read and mapping statistics (from step 3)
- Nextclade clade/lineage assignments and genome coverage (from step 5)
- UShER closest-neighbour results: distance, number of tied neighbours, countries, and date range (from step 12)
- Sample metadata from the samplesheet

A sequence quality classification is assigned based on genome coverage and lineage assignment status:

| Coverage | Lineage | Quality text |
|---|---|---|
| >80% | Assigned | `High Quality. Lineage and phylogenetic placement high confidence.` |
| >80% | Unassigned/NA | `High Quality genome. Lineage unassigned. Phylogenetic placement high confidence.` |
| 50–80% | Assigned | `Medium Quality. Lineage and phylogenetic placement medium confidence.` |
| 50–80% | Unassigned/NA | `Medium Quality genome. Lineage unassigned. Phylogenetic placement medium confidence.` |
| <50% | Assigned | `Low Quality. Lineage and phylogenetic placement low confidence.` |
| <50% | Unassigned/NA | `Low Quality genome. Lineage unassigned. Phylogenetic placement low confidence.` |

**Input:** Combined read stats TSV, combined depth stats TSV, samplesheet, Nextclade CSV, closest-neighbour report TSV  
**Output:** `{RunName}_final_results.csv`

---

## Output folder structure

All results are written to `params.output_dir/{RunName}/`:

```
{RunName}/
├── filtering/
│   └── {barcode}_{sample_id}_filtered.fastq       # Length-filtered reads per sample
│
├── minion/
│   ├── {barcode}_{sample_id}/                      # All artic minion output files per sample
│   │   ├── {barcode}_{sample_id}.consensus.fasta   # Consensus genome
│   │   ├── {barcode}_{sample_id}.primertrimmed.rg.sorted.bam  # Primer-trimmed BAM
│   │   ├── {barcode}_{sample_id}.sorted.bam        # Un-normalised BAM
│   │   └── ...                                     # Other artic minion intermediates
│   ├── read_stats/
│   │   ├── {barcode}_{sample_id}_read_stats.tsv    # Per-sample read/mapping counts
│   │   └── {RunName}_combined_read_stats.tsv       # All samples combined
│   └── depth_stats/
│       ├── {barcode}_{sample_id}_depth_stats.tsv   # Per-sample depth statistics
│       └── {RunName}_combined_depth_stats.tsv      # All samples combined
│
├── combined/
│   └── {RunName}_combined_consensus.fasta          # All consensus genomes concatenated
│
├── nextclade_output/
│   ├── {RunName}.csv                               # Per-sample clade, lineage, coverage
│   └── ...                                         # Full Nextclade output folder
│
├── {RunName}_final_results.csv                     # Master summary (see below)
│
└── usher/
    ├── {RunName}_full_tree_global.nwk              # Newick: global placement tree
    ├── {RunName}_full_tree_global.json             # JSON: global tree with global metadata
    ├── {RunName}_full_tree_optimized.nwk           # Newick: optimised tree
    ├── {RunName}_full_tree_optimized.json          # JSON: optimised tree with global metadata
    ├── {RunName}_context_tree.nwk                  # Newick: context subtree (20 neighbours)
    ├── {RunName}_context_tree.json                 # JSON: context subtree with metadata
    ├── {RunName}_final_report.tsv                  # Closest global neighbour per sample
    ├── {RunName}_final_report_10SNPs.tsv           # All neighbours within ≤10 SNPs
    └── {RunName}_closest_neighbor_report.tsv       # R summary report (see below)
```

---

## Key output files

### `{RunName}_final_results.csv`
Master per-sample summary table produced by the R summary script. Contains:

| Column | Description |
|---|---|
| `LW_id` | Laboratory sample ID |
| `Barcode` | Sequencing barcode |
| `Clade` | MPXV clade (Nextclade) |
| `Lineage` | Pango lineage (Nextclade); blank if unassigned |
| `Coverage` | Genome coverage (%) |
| `RawReads` | Total reads before filtering |
| `TrimmedReads` | Reads after length filtering |
| `MappedReads` | Reads mapped to reference (primer-trimmed BAM) |
| `MappedPercentage` | Mapping rate (%) |
| `AverageDepth` | Mean sequencing depth (un-normalised BAM) |
| `Phylo_Distance` | Parsimony distance to closest global neighbour(s) |
| `N_Tied_Neighbors` | Number of equally close global neighbours |
| `Neighbor_Countries` | Countries of closest global neighbour(s) |
| `Neighbor_Date_Range` | Collection date range of closest global neighbour(s) |
| `Sequence quality` | Quality classification with lineage and phylogenetic placement confidence (see Step 6) |

---

## UShER Output Files (Most Relevant)

The following outputs from the UShER phylogenetic placement workflow are most useful for downstream analysis:

### `{RunName}_full_tree_optimized.nwk`
Optimised phylogenetic tree in Newick format. Contains all global MPXV sequences plus the run samples with optimised topology. Recommended for publication-quality phylogenetic figures using R (`ggtree`, `ape`) or other phylogenetic tools.

---

### `{RunName}_context_tree.nwk` and `{RunName}_context_tree.json`
Context subtree retaining the 20 closest global neighbours per local sample. The `.nwk` version is ideal for focused visualisation and figure generation. The `.json` version is compatible with [Taxonium](https://taxonium.org/) or [Auspice](https://auspice.us/) for interactive phylogenetic exploration and visualization of closely related sequences.

---

### `{RunName}_closest_neighbor_report.tsv`
Tab-separated summary table with one row per local sample:

| Column | Description |
|---|---|
| `My_Sample` | Local sample tip label (`{barcode}_{sample_id}`) |
| `Count_Neighbors` | Number of equally close global neighbours (ties) |
| `Distance` | Parsimony distance to closest neighbour(s) |
| `All_Closest_Neighbors` | Semicolon-separated list of closest global strain names |
| `Countries` | Pipe-separated list of countries of closest neighbours |
| `Date_Range` | Date range of closest neighbours (`YYYY-MM-DD` or range) |

---

## Additional UShER Output Files

### `{RunName}_final_report.tsv`
Per-sample closest global neighbour report produced by `matUtils`. Lists the single most parsimonious placement relative in the global tree for each local sample, including parsimony distance and strain name.

---

### `{RunName}_final_report_10SNPs.tsv`
Extended nearest-neighbour report listing all global sequences within ≤10 SNPs (parsimony distance) of each local sample. Useful for outbreak cluster investigation.

---

### `{RunName}_full_tree_optimized.json` and `{RunName}_full_tree_global.json`
JSON tree files compatible with [Taxonium](https://taxonium.org/) or [Auspice](https://auspice.us/) for interactive phylogenetic visualisation. The full tree contains all global sequences plus the run samples. Use the optimised version for best topology.

