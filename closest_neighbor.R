library(tidyverse)
library(ggtree)
library(treeio)
library(tidytree)
library(ape)
library(castor)


#Import global meta
global_meta <- read_tsv("../UShER/mpxv.cladeII.2026-01-09.metadata.tsv") %>% 
  mutate(sample_type = "Global",
         date = as.character(date)) %>% 
  select(strain, date, country, sample_type)

#Because one or more of the global samples contain a ' in the sample name we need to remove this from the name before loading in the treefiles....

# 1. Read the file as raw text first
raw_text <- readLines("../82_samples_UShER_matoptimize_v7/full_tree_optimized.nwk", warn = FALSE)

# 2. Remove the offending single quote (replace ' with nothing)
#    This turns "...Belgium/2022'|OR..." into "...Belgium/2022|OR..."
clean_text <- gsub("'", "", raw_text)

# 3. Now parse the cleaned text into a tree object
phylo_global <- read.tree(text = clean_text)

#Import own metadata
#THIS NEEDS TO CHANGE. 
#IT NEEDS A PART THAT TAKES THE HEADER FROM THE FASTA FILES AS SAMPLENAMES AND THE SAMPLEDATES FROM THE SAMPLESHEET FROM THE CORRECT SAMPLE AND CREATES THE SAMPLES META FROM SCRATCH
#HEADER EXAMPLE: >barcode65_1037456 NC_063383_artic-network/fieldbioinformatics_1.8.5 
#Tip.lab example: barcode65_1037456 NC_063383_artic-network/fieldbioinformatics_1.8.5 #WOULD THIS KEEP THE PART AFTE THE SPACE? UNSURE. 
samples_meta <- read_tsv("../filtered_sequences/filtered_samples_headername_dates.tsv") %>% 
  rename(strain = sample) %>% 
  mutate(country = "Norway",
         sample_type = "Project",
         date = as.character(date)) 
#Combine global and local metadata
combined_meta <- bind_rows(samples_meta, global_meta) %>% 
  rename(label = strain)

#Add metadata to treefile
tree_with_data <- full_join(phylo_global, combined_meta, by = "label")

#Find closes neighbor and create final report table 
# 1. Get the names and Match indices (Same as before)
global_samples <- combined_meta %>% filter(sample_type != "Project") %>% pull(label)
my_samples <- combined_meta %>% filter(sample_type == "Project") %>% pull(label)

global_indices <- match(global_samples, phylo_global$tip.label)
my_indices <- match(my_samples, phylo_global$tip.label)

# Filter valid indices
valid_global_indices <- global_indices[!is.na(global_indices)]
clean_my_indices <- my_indices[!is.na(my_indices)]
clean_my_labels <- my_samples[!is.na(my_indices)]

# 2. Setup Loop
results_list <- list()
tolerance <- 1e-9 # Tolerance for "equal" distance

message(paste("Processing", length(clean_my_indices), "samples..."))

for (i in seq_along(clean_my_indices)) {
  
  target_idx <- clean_my_indices[i]
  target_name <- clean_my_labels[i]
  
  # --- FIX USING YOUR SCREENSHOT FUNCTION ---
  # Calculate distance from THIS sample to ALL other tips
  # Returns a vector where index 1 = Tip 1, index 2 = Tip 2, etc.
  all_dists <- get_all_distances_to_tip(phylo_global, focal_tip = target_idx)
  
  # ------------------------------------------
  
  # B. Subset to just the Global samples
  # We extract the distances corresponding to the global sample indices
  global_dists <- all_dists[valid_global_indices]
  
  # C. Find the minimum distance
  min_dist <- min(global_dists)
  
  # D. Find ALL global indices that match this minimum
  is_closest <- abs(global_dists - min_dist) < tolerance
  
  tie_indices <- valid_global_indices[is_closest]
  
  # E. Store results
  results_list[[i]] <- data.frame(
    My_Sample = target_name,
    Closest_Global_Sample = phylo_global$tip.label[tie_indices],
    Distance = min_dist
  )
}


# 3. Combine results
final_report_ties <- bind_rows(results_list)

# 4. Join with metadata
# Note: This will now work because 'final_report_ties' is properly populated
final_report_ties <- final_report_ties %>%
  left_join(combined_meta, by = c("Closest_Global_Sample" = "label")) %>%
  select(My_Sample, Closest_Global_Sample, Distance, country, date) %>%
  arrange(My_Sample)

# 5. View Result
print(head(final_report_ties))


library(dplyr)
library(lubridate)

summary_report <- final_report_ties %>%
  # Fix dates (allow partial dates like "2022-05")
  mutate(date_clean = ymd(date, quiet = TRUE)) %>%
  
  group_by(My_Sample) %>%
  
  summarise(
    # 1. NEW COLUMN: Count the number of tied neighbors
    Count_Neighbors = n(),
    
    # 2. Distance (same for all)
    Distance = first(Distance),
    
    # 3. List names
    All_Closest_Neighbors = paste(Closest_Global_Sample, collapse = "; "),
    
    # 4. List Countries
    Countries = paste(unique(na.omit(country)), collapse = "; "),
    
    # 5. Date Range
    Min_Date = as.character(min(date_clean, na.rm = TRUE)),
    Max_Date = as.character(max(date_clean, na.rm = TRUE))
  ) %>%
  
  # Format Range
  mutate(Date_Range = ifelse(Min_Date == Max_Date, 
                             Min_Date, 
                             paste(Min_Date, "to", Max_Date))) %>%
  
  ungroup()

# View the result
head(summary_report)
