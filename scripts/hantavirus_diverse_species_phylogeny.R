library(ggtree)
library(ape)
library(tidyverse)
library(phangorn)

tree <- read.tree("data/global_hantavirus_species_phylogeny/phylogeny.treefile") %>% 
  midpoint
metadata <- read.csv("data/global_hantavirus_species_phylogeny/orthohantaviruses_ncbi_virus_metadata.csv")
metadata <- metadata %>%
  mutate(Accession = paste0(Accession, ".1"))
metadata <- metadata %>% 
  mutate(is_outbreak = if_else(Accession == "PZ385162.1", "Outbreak", "Reference"),
         Geo_Location == if_else(Geo_Location == "" | is.na(Geo_Location), "NA",  Geo_Location))

andes_node <- getMRCA(tree, tip = c("PZ385162.1", "PV808474.1", "MN258191.1"))

p <- ggtree(tree) %<+% metadata +
  
  geom_tiplab(aes(label = Species), size = 4, offset = 0.005)+
  
  geom_tippoint(data = ~ filter(.x, is_outbreak == "Outbreak"),
                color = "red", shape = 18, size = 7)+
  
  geom_tippoint(aes(color = Geo_Location), size = 2.5, alpha = 0.8)+
  
  geom_hilight(node = andes_node, fill = "steelblue", alpha = 0.2, extend = 1.5) +  
  
  hexpand(0.5)
  
p


library(ggtree)
library(ape)
library(tidyverse)
library(phangorn)

# 1. Load Data with Relative Paths
tree <- read.tree("data/global_hantavirus_species_phylogeny/phylogeny.treefile") %>% 
  midpoint()

metadata <- read.csv("data/global_hantavirus_species_phylogeny/orthohantaviruses_ncbi_virus_metadata.csv")

# 2. Clean Metadata
metadata <- metadata %>%
  mutate(Accession = paste0(Accession, ".1")) %>%
  mutate(
    is_outbreak = if_else(Accession == "PZ385162.1", "Outbreak", "Reference"),
    # FIXED: Changed == to = and removed the tilde (~)
    Geo_Location = if_else(Geo_Location == "" | is.na(Geo_Location), "NA", Geo_Location)
  )

# 3. Find the Andes Node
# Ensure these accessions exist in your tree tips exactly!
#andes_node <- getMRCA(tree, tip = c("PZ385162.1", "PV808474.1", "MN258191.1"))
outbreak_node <- which(tree$tip.label == "PZ385162.1")


# 4. Generate Plot
p <- ggtree(tree) %<+% metadata +
  geom_tiplab(aes(label = Species), size = 3, offset = 0.5, align = TRUE) +
  geom_tippoint(data = ~ filter(.x, is_outbreak == "Outbreak"),
                color = "red", shape = 18, size = 6) +
  geom_tippoint(aes(color = Geo_Location), size = 2.5, alpha = 0.8) +
  geom_hilight(node = outbreak_node, fill = "firebrick", alpha = 0.2, extend = 4.2) +  
  geom_text2(aes(label = label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 95), 
             vjust = -0.5, 
             hjust = 1.1, 
             size = 3, 
             color = "grey30")+
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)))+
  guides(color = guide_legend(ncol = 3, title = "Country", override.aes = list(size = 4)))+
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9,  face = "bold"),
        legend.box.margin = margin(t=2),
        plot.margin = margin(10,10,10,10))+
  hexpand(0.1)

# 5. Return the plot object
p