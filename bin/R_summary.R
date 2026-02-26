library(tidyverse)
options <- commandArgs(trailingOnly = TRUE)

read_stats   <- read_tsv(options[1])
sample_sheet <- read_csv2(options[2])
nextclade    <- read_csv2(options[3])
depth_stats  <- read_tsv(options[5])
neighbor_report <- read_tsv(options[6])

sample_sheet <- sample_sheet %>%
  mutate(
    sample_id = as.character(sample_id),
    barcode   = as.character(barcode)
  )


read_stats <- read_stats %>% 
  rename(sample_id = Sample) %>%
  separate(sample_id, into = c("barcode", "sample_id"), sep = "_") %>%
  mutate(
    sample_id = as.character(sample_id),
    barcode   = as.character(barcode)
  )

#Clean up the seqName column to only retain barcode and sample ID and massage the data
nextclade$seqName <- sub(" .*", "", nextclade$seqName)
nextclade <- nextclade %>%
  rename(sample_id = seqName) %>%
  select(sample_id, clade, lineage, coverage) %>%
  mutate(
    coverage = as.numeric(coverage) * 100
  ) %>%
  separate(sample_id, into = c("barcode", "sample_id"), sep = "_") %>%
  mutate(
    sample_id = as.character(sample_id),
    barcode   = as.character(barcode)
  )

depth_stats <- depth_stats %>%
  rename(sample_id = Sample) %>%
  separate(sample_id, into = c("barcode", "sample_id"), sep = "_") %>%
  mutate(
    sample_id = as.character(sample_id),
    barcode   = as.character(barcode)
  ) %>%
  # keep only the three new columns to avoid name clashes
  select(
    barcode, sample_id,
    Mapped_Reads_Unnormalized,
    Mapped_Percentage_Unnormalized,
    AverageDepth
  )

combined_all <- sample_sheet %>%
  left_join(nextclade,  by = c("sample_id", "barcode")) %>%
  left_join(read_stats, by = c("sample_id", "barcode")) %>%
  left_join(depth_stats, by = c("sample_id", "barcode"))

# Join UShER closest-neighbor columns using reconstructed tip label
neighbor_report <- neighbor_report %>%
  rename(tip_label = My_Sample) %>%
  select(tip_label, Distance, Count_Neighbors, Countries, Date_Range)

combined_all <- combined_all %>%
  mutate(tip_label = paste0(barcode, "_", sample_id)) %>%
  left_join(neighbor_report, by = "tip_label") %>%
  select(-tip_label)

combined_all$Species <- "M-koppevirus"

combined_all <- combined_all %>% 
  rename(
    LW_id = sample_id,
    Barcode = barcode, 
    Clade = clade, 
    Lineage = lineage, 
    Coverage = coverage, 
    RawReads = Raw_Reads,
    TrimmedReads = Filtered_Reads, 
    MappedReads = Mapped_Reads, 
    MappedPercentage = Mapped_Percentage,
    Phylo_Distance = Distance,
    N_Tied_Neighbors = Count_Neighbors,
    Neighbor_Countries = Countries,
    Neighbor_Date_Range = Date_Range
  ) %>% 
  select(
    -fastq, -Mapped_Reads_Unnormalized, -Mapped_Percentage_Unnormalized
  ) %>%
  # Replace semicolon separators in Neighbor_Countries with pipe separators for proper CSV formatting
  mutate(Neighbor_Countries = gsub("; ", " | ", Neighbor_Countries, fixed = TRUE)) %>%
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
