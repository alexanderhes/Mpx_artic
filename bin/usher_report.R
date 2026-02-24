library(tidyverse)
library(lubridate)
library(ape)
library(castor)

# ==============================================================================
# Arguments
#   1: combined consensus FASTA (pipeline output from COMBINE_FASTA)
#   2: filled samplesheet (semicolon-separated, contains SampleDate column)
#   3: optimized NWK tree (pipeline output from USHER_PLACE)
#   4: global metadata TSV (downloaded by USHER_DOWNLOAD)
#   5: output CSV path
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: usher_report.R <fasta> <samplesheet> <tree.nwk> <global_metadata.tsv> <output.csv>")
}

fasta_file       <- args[1]
samplesheet_file <- args[2]
tree_file        <- args[3]
global_meta_file <- args[4]
output_file      <- args[5]

# ==============================================================================
# 1. Extract tip labels from FASTA headers
#    Header format: >barcode65_1037456 NC_063383_artic-network/fieldbioinformatics_1.8.5
#    Tip label:      barcode65_1037456  (first word only, no ">")
# ==============================================================================
fasta_lines  <- readLines(fasta_file)
header_lines <- fasta_lines[startsWith(fasta_lines, ">")]
tip_labels   <- sub("^>([^ ]+).*", "\\1", header_lines)

message(paste("Found", length(tip_labels), "sequences in FASTA"))

# ==============================================================================
# 2. Build local sample metadata from samplesheet
#    - Construct strain ID as {barcode}_{sample_id} to match FASTA headers
#    - Parse SampleDate (Norwegian DD.MM.YYYY format)
# ==============================================================================
samplesheet <- read_csv2(samplesheet_file, show_col_types = FALSE) %>%
  mutate(
    sample_id = as.character(sample_id),
    barcode   = tolower(as.character(barcode)),
    strain    = paste0(barcode, "_", sample_id),
    date      = format(dmy(SampleDate), "%Y-%m-%d"),
    country   = "Norway",
    sample_type = "Project"
  ) %>%
  filter(strain %in% tip_labels) %>%
  select(strain, date, country, sample_type)

if (nrow(samplesheet) == 0) {
  stop("No samplesheet entries matched FASTA tip labels. Check barcode/sample_id formatting.")
}

message(paste("Matched", nrow(samplesheet), "local samples from samplesheet"))

# ==============================================================================
# 3. Load global metadata
# ==============================================================================
global_meta <- read_tsv(global_meta_file, show_col_types = FALSE) %>%
  mutate(
    sample_type = "Global",
    date        = as.character(date)
  ) %>%
  select(strain, date, country, sample_type)

# ==============================================================================
# 4. Combine local and global metadata
# ==============================================================================
combined_meta <- bind_rows(samplesheet, global_meta) %>%
  rename(label = strain)

# ==============================================================================
# 5. Load and clean tree
#    Some global sample names contain a single quote which breaks read.tree()
# ==============================================================================
raw_text   <- readLines(tree_file, warn = FALSE)
clean_text <- gsub("'", "", raw_text)
phylo_tree <- read.tree(text = clean_text)

message(paste("Tree loaded:", ape::Ntip(phylo_tree), "tips"))

# ==============================================================================
# 6. Find closest global neighbor for each local sample
# ==============================================================================
global_samples <- combined_meta %>% filter(sample_type == "Global") %>% pull(label)
my_samples     <- combined_meta %>% filter(sample_type == "Project") %>% pull(label)

global_indices <- match(global_samples, phylo_tree$tip.label)
my_indices     <- match(my_samples,     phylo_tree$tip.label)

valid_global_indices <- global_indices[!is.na(global_indices)]
clean_my_indices     <- my_indices[!is.na(my_indices)]
clean_my_labels      <- my_samples[!is.na(my_indices)]

tolerance    <- 1e-9
results_list <- list()

message(paste("Processing", length(clean_my_indices), "local samples..."))

for (i in seq_along(clean_my_indices)) {

  target_idx  <- clean_my_indices[i]
  target_name <- clean_my_labels[i]

  # Distance from this tip to all other tips in the tree
  all_dists    <- get_all_distances_to_tip(phylo_tree, focal_tip = target_idx)
  global_dists <- all_dists[valid_global_indices]

  min_dist   <- min(global_dists)
  is_closest <- abs(global_dists - min_dist) < tolerance
  tie_indices <- valid_global_indices[is_closest]

  results_list[[i]] <- data.frame(
    My_Sample            = target_name,
    Closest_Global_Sample = phylo_tree$tip.label[tie_indices],
    Distance             = min_dist
  )
}

final_report <- bind_rows(results_list) %>%
  left_join(combined_meta, by = c("Closest_Global_Sample" = "label")) %>%
  select(My_Sample, Closest_Global_Sample, Distance, country, date) %>%
  arrange(My_Sample)

# ==============================================================================
# 7. Summarise: one row per local sample
# ==============================================================================
summary_report <- final_report %>%
  mutate(date_clean = ymd(date, quiet = TRUE)) %>%
  group_by(My_Sample) %>%
  summarise(
    Count_Neighbors       = n(),
    Distance              = first(Distance),
    All_Closest_Neighbors = paste(Closest_Global_Sample, collapse = "; "),
    Countries             = paste(unique(na.omit(country)), collapse = "; "),
    Min_Date              = as.character(min(date_clean, na.rm = TRUE)),
    Max_Date              = as.character(max(date_clean, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    Date_Range = ifelse(Min_Date == Max_Date, Min_Date, paste(Min_Date, "to", Max_Date))
  ) %>%
  select(-Min_Date, -Max_Date)

message(paste("Writing report to:", output_file))
write_csv(summary_report, output_file)
message("Done.")
