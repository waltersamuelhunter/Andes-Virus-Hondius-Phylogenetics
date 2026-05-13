library(tidytree)
library(treeio)
library(patchwork)
library(ggtree)
library(cowplot)


cruise_isolates <- c("PP_006WDKH.1", "PP_006WDJK.1", "PP_006WANE.2", "PP_006W6RC.2", "PP_006WBLH.1", "PP_006W3U9.2")

cruise_outbreak <- getMRCA(tree_M, tip = cruise_isolates)

metadata_M <- read.csv("M_segment_pathoplexus_andes_26_05_13/andv_141seq_pathoplexus_meta.tsv" , sep = "\t")
tree_M <- read.tree("M_segment_pathoplexus_andes_26_05_13/tree_result/phylogeny.treefile") %>% midpoint 
mrca_M <- getMRCA(tree_M, tip = c("PP_006WDKH.1", "PP_006WDJK.1", "PP_006WANE.2", "PP_006W6RC.2", "PP_006WBLH.1", "PP_006W3U9.2"))
p_node_M <- parent(as_tibble(tree_M), mrca_M)$node
gp_node_M <- parent(as_tibble(tree_M), p_node_M)$node
M_subtree <- extract.clade(tree_M, node = gp_node_M)

metadata_S <- read.csv("S_segment_pathoplexus_andes_26_05_13/andv_metadata_2026-05-13T0723.tsv", sep = "\t")
tree_S <- read.tree("S_segment_pathoplexus_andes_26_05_13/tree_results/phylogeny.treefile") %>% midpoint
mrca_S <- getMRCA(tree_S, tip = c("PP_006WDKH.1", "PP_006WDJK.1", "PP_006WANE.2", "PP_006W6RC.2", "PP_006WBLH.1", "PP_006W3U9.2"))
p_node_S <- parent(as_tibble(tree_S), mrca_S)$node
gp_node_S <- parent(as_tibble(tree_S), p_node_S)$node
S_subtree <- extract.clade(tree_S, node = gp_node_S)

metadata_L <- read.csv("L_segment_pathoplexus_andes_26_05_13/andv_metadata_2026-05-13T0726.tsv", sep = "\t")
tree_L <- read.tree("L_segment_pathoplexus_andes_26_05_13/tree_results/phylogeny.treefile") %>% midpoint
mrca_L <- getMRCA(tree_L, tip = c("PP_006WDKH.1", "PP_006WDJK.1", "PP_006WANE.2", "PP_006W6RC.2", "PP_006WBLH.1", "PP_006W3U9.2"))
p_node_L <- parent(as_tibble(tree_L), mrca_L)$node
gp_node_L <- parent(as_tibble(tree_L), p_node_L)$node
L_subtree <- extract.clade(tree_L, node = gp_node_L)


all_countries <- sort(unique(c(metadata_S$geoLocCountry, 
                               metadata_M$geoLocCountry, 
                               metadata_L$geoLocCountry)))

plot_subtree <- function(subtree, metadata, highlight_isolates, title_text) {

  target_node <- getMRCA(subtree, tip = intersect(highlight_isolates, subtree$tip.label))
  
  p <-  ggtree(subtree) %<+% metadata + 
  
  geom_hilight(node = target_node,
                 fill = "firebrick",
                 extend = 0.11,
                 alpha = 0.2,
                 show.legend = FALSE) +
  
  geom_tippoint(aes(fill = geoLocCountry,
                    shape = hostNameScientific),
                size = 3.5, 
                color = "black",
                show.legend = FALSE)+
  
  geom_tiplab(size = 2.5,
              offset = 0.03)+
  
  scale_fill_viridis_d(option = "plasma", 
                       na.value = "grey90", 
                       name = "Country",
                       limits = all_countries,
                       na.translate = FALSE) +
    
  scale_shape_manual(values = c(21,24,22,23,25,8), name = "Host") +
 
  theme_tree2() +

  labs(title = title_text)+
  
  hexpand(2.5)+
  
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 5), order = 1),
           shape = guide_legend(override.aes = list(size = 4), order = 2))
  
}


p1 <- plot_subtree(subtree = S_subtree, 
                   metadata = metadata_S, 
                   highlight_isolates = cruise_isolates,
                   title_text = "segment_S")

p2 <- plot_subtree(subtree = M_subtree, 
                   metadata = metadata_M, 
                   highlight_isolates = cruise_isolates,
                   title_text = "segment_M")

p3 <- plot_subtree(subtree = L_subtree, 
                   metadata = metadata_L, 
                   highlight_isolates = cruise_isolates,
                   title_text = "segment_L")


# Create a dummy plot to generate the legend
legend_plot <- ggplot(metadata_M, aes(x=1, y=1, fill=geoLocCountry, shape=hostNameScientific)) +
  geom_point(size = 4, color = "black", stroke = 0.4) +
  scale_fill_viridis_d(option = "plasma", name = "Country", limits = all_countries) +
  scale_shape_manual(values = c(21, 24, 22, 23, 25, 8), name = "Host Species") +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 5)),
         shape = guide_legend(override.aes = list(size = 4))) +
  theme_void() +
  theme(legend.position = "right")

# Extract the legend using cowplot or patchwork

shared_legend <- get_legend(legend_plot)

final_plot <- (p1 | p2 | p3) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.title = element_text(size = 9, face = "bold"),
        plot.title = element_text(size = 10, face = "bold"))

# p1, p2, p3 now have show.legend = FALSE
combined_trees <- (p1 | p2 | p3) 

# Final assembly
final_output <- plot_grid(combined_trees, shared_legend, rel_widths = c(3, 0.6))

final_output

