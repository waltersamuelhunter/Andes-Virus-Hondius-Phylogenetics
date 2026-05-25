# scripts/extract_tmrca_function.R

library(ape)
library(ggtree)
library(phangorn)
library(lubridate)

# Define the custom function
extract_outbreak_dates <- function(nexus_path, isolates, mrsd = "2026-05-06") {
  
  # 1. Load and root the tree
  time_tree <- read.nexus(nexus_path) %>% midpoint()
  
  # 2. Find the Crown node (MRCA)
  crown_node <- getMRCA(time_tree, isolates)
  
  # 3. Build the tree data object
  p <- ggtree(time_tree, mrsd = mrsd)
  tree_df <- p$data
  
  # 4. Find the Stem node (the immediate parent of the Crown node)
  stem_node <- tree_df$parent[tree_df$node == crown_node]
  
  # 5. Extract and convert the dates
  crown_clean <- as.numeric(tree_df[tree_df$node == crown_node, "x"])
  crown_date <- as.character(as.Date(date_decimal(crown_clean)))
  
  stem_clean <- as.numeric(tree_df[tree_df$node == stem_node, "x"])
  stem_date <- as.character(as.Date(date_decimal(stem_clean)))
  
  # Return both dates as a list
  return(list(Stem = stem_date, Crown = crown_date))
}