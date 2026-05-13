library(phangorn)
library(seqinr)

metadata_M <- read.csv("M_segment_pathoplexus_andes_26_05_13/andv_141seq_pathoplexus_meta.tsv" , sep = "\t")
tree_M <- read.tree("M_segment_pathoplexus_andes_26_05_13/tree_result/phylogeny.treefile") %>% midpoint 
aln_M <- read.fasta("M_segment_pathoplexus_andes_26_05_13/aln.fasta")

metadata_S <- read.csv("S_segment_pathoplexus_andes_26_05_13/andv_metadata_2026-05-13T0723.tsv", sep = "\t")
tree_S <- read.tree("S_segment_pathoplexus_andes_26_05_13/tree_results/phylogeny.treefile") %>% midpoint
aln_S <- read.fasta("S_segment_pathoplexus_andes_26_05_13/andv_aligned-nuc-S_2026-05-13T0724.fasta")

metadata_L <- read.csv("L_segment_pathoplexus_andes_26_05_13/andv_metadata_2026-05-13T0726.tsv", sep = "\t")
tree_L <- read.tree("L_segment_pathoplexus_andes_26_05_13/tree_results/phylogeny.treefile") %>% midpoint
aln_L <- read.fasta("L_segment_pathoplexus_andes_26_05_13/andv_aligned-nuc-L_2026-05-13T0726.fasta")


find_stem_snps <- function(aln, outbreak_isolates, sister_group) {
  # 1. Subset the alignment list
  out_list <- aln[names(aln) %in% outbreak_isolates]
  sis_list <- aln[names(aln) %in% sister_group]
  
  # 2. Convert list to character matrix (each row is a sequence)
  # We use do.call(rbind, ...) to stack the character vectors
  out_mat <- do.call(rbind, out_list)
  sis_mat <- do.call(rbind, sis_list)
  
  # 3. Logic: Find positions where:
  # - All outbreak seqs are identical to each other
  # - All sister seqs are identical to each other
  # - Outbreak differs from Sister
  
  # Identify positions (columns)
  is_fixed_outbreak <- apply(out_mat, 2, function(x) length(unique(x)) == 1)
  is_fixed_sister   <- apply(sis_mat, 2, function(x) length(unique(x)) == 1)
  is_different      <- out_mat[1, ] != sis_mat[1, ]
  
  snp_indices <- which(is_fixed_outbreak & is_fixed_sister & is_different)
  
  if(length(snp_indices) == 0) return(list(positions = NULL, labels = "No fixed SNPs"))
  
  # 4. Create labels (e.g., "A124G")
  labels <- sapply(snp_indices, function(i) {
    paste0(toupper(sis_mat[1, i]), i, toupper(out_mat[1, i]))
  })
  
  return(list(positions = snp_indices, labels = paste(labels, collapse = ", ")))
}

# Run for all three segments
snps_M <- find_stem_snps(aln_M, cruise_isolates, sister)
snps_S <- find_stem_snps(aln_S, cruise_isolates, sister)
snps_L <- find_stem_snps(aln_L, cruise_isolates, sister)
cruise_isolates <- c("PP_006WDKH.1", "PP_006WDJK.1", "PP_006WANE.2", "PP_006W6RC.2", "PP_006WBLH.1", "PP_006W3U9.2")
sister <- c("PP_006VZ5U.1", "PP_006VZ4W.1")


snps_M <- find_stem_snps(aln_M, cruise_isolates, sister)
snps_S <- find_stem_snps(aln_S, cruise_isolates, sister)
snps_L <- find_stem_snps(aln_L, cruise_isolates, sister)

# Translate a specific SNP to check for synonymous vs non-synonymous
# This assumes your alignment is in-frame
check_mutation_type <- function(aln_mat, pos) {
  # Logic to find codon, translate, and compare
  # Useful for your GitHub 'Technical Deep Dive' section
}

check_mutation_type <- function(aln, pos, outbreak_isolates, sister_group) {
  # 1. Identify the codon position (assuming alignment starts at 1st codon)
  codon_start <- ((pos - 1) %/% 3) * 3 + 1
  codon_indices <- codon_start:(codon_start + 2)
  
  # 2. Extract sequences as matrix
  mat <- do.call(rbind, aln)
  
  # 3. Get representative codons
  out_codon <- toupper(paste(mat[which(names(aln) %in% outbreak_isolates)[1], codon_indices], collapse = ""))
  sis_codon <- toupper(paste(mat[which(names(aln) %in% sister_group)[1], codon_indices], collapse = ""))
  
  # 4. Translate using seqinr
  out_aa <- translate(s2c(out_codon))
  sis_aa <- translate(s2c(sis_codon))
  
  # 5. Result
  type <- if(out_aa == sis_aa) "Synonymous (S)" else "Non-synonymous (NS)"
  
  return(data.frame(
    Pos = pos,
    Codon_Pos = (pos - 1) %/% 3 + 1,
    Sister = sis_aa,
    Outbreak = out_aa,
    Type = type,
    label = paste0(sis_aa, ((pos - 1) %/% 3 + 1), out_aa)
  ))
}

# Example for M segment
m_mutations <- lapply(snps_M$positions, function(p) {
  check_mutation_type(aln_M, p, cruise_isolates, sister)
})

m_summary <- do.call(rbind, m_mutations)
print(m_summary)