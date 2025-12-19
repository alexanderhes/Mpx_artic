library(tidyverse)
options <- commandArgs(trailingOnly = TRUE)

read_stats <- read_tsv(options[1])
sample_sheet <- read_csv2(options[2])
nextclade <- read_csv2(options[3])
depth_stats  <- read_tsv(options[5])

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

combined_all$Species <- "M-koppevirus"

combined_all <- combined_all %>% 
  rename(
    SequenceID = sample_id,
    Barcode = barcode, 
    Clade = clade, 
    Lineage = lineage, 
    Coverage = coverage, 
    RawReads = Raw_Reads,
    TrimmedReads = Filtered_Reads, 
    MappedReads = Mapped_Reads, 
    MappedPercentage = Mapped_Percentage
  ) %>% 
  select(
    -fastq, -PrøveNr
  )

write_csv(combined_all, options[4])
