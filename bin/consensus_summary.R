library(tidyverse)

# ==============================================================================
# Final results for consensus mode (no read/depth stats)
#   1: Nextclade CSV
#   2: prepared samplesheet (from PREPARE_CONSENSUS, contains tip_label)
#   3: closest neighbor report TSV (from USHER_REPORT)
#   4: output CSV path
# ==============================================================================
options <- commandArgs(trailingOnly = TRUE)

# Read everything as text so sample IDs keep leading zeros
as_text <- cols(.default = col_character())
nextclade       <- read_csv2(options[1], col_types = as_text)
sample_sheet    <- read_csv2(options[2], col_types = as_text)
neighbor_report <- read_tsv(options[3],  col_types = as_text)

nextclade <- nextclade %>%
  transmute(
    tip_label = sub(" .*", "", seqName),
    Clade     = clade,
    Lineage   = lineage,
    Coverage  = as.numeric(coverage) * 100,
    QC_Status = `qc.overallStatus`
  )

neighbor_report <- neighbor_report %>%
  transmute(
    tip_label           = My_Sample,
    Phylo_Distance      = Distance,
    N_Tied_Neighbors    = Count_Neighbors,
    # Pipe separators for proper CSV formatting
    Neighbor_Countries  = gsub("; ", " | ", Countries, fixed = TRUE),
    Neighbor_Date_Range = Date_Range
  )

combined_all <- sample_sheet %>%
  left_join(nextclade,       by = "tip_label") %>%
  left_join(neighbor_report, by = "tip_label") %>%
  select(-tip_label) %>%
  rename(LW_id = sample_id) %>%
  mutate(Species = "M-koppevirus") %>%
  # Sequence quality: coverage-based confidence + lineage assignment status
  mutate(`Sequence quality` = case_when(
    Coverage > 80  & !is.na(Lineage) & Lineage != "unassigned" ~ "High Quality. Lineage and phylogenetic placement high confidence.",
    Coverage > 80  & (is.na(Lineage) | Lineage == "unassigned") ~ "High Quality genome. Lineage unassigned. Phylogenetic placement high confidence.",
    Coverage >= 50 & Coverage <= 80 & !is.na(Lineage) & Lineage != "unassigned" ~ "Medium Quality. Lineage and phylogenetic placement medium confidence.",
    Coverage >= 50 & Coverage <= 80 & (is.na(Lineage) | Lineage == "unassigned") ~ "Medium Quality genome. Lineage unassigned. Phylogenetic placement medium confidence.",
    Coverage < 50  & !is.na(Lineage) & Lineage != "unassigned" ~ "Low Quality. Lineage and phylogenetic placement low confidence.",
    Coverage < 50  & (is.na(Lineage) | Lineage == "unassigned") ~ "Low Quality genome. Lineage unassigned. Phylogenetic placement low confidence.",
    TRUE           ~ "Unknown"
  ))

write_csv(combined_all, options[4])
