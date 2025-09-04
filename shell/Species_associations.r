# BioDiv_Vienna.r
# Analysis of Drosophila biodiversity in Vienna
# Author: [Your Name]
# Date: [Date]
# Description: This script performs PCA, diversity estimation, NMDS ordination, and mixed-effects modeling on Vienna Drosophila samples.

# --- Load Required Libraries ---

# --- Set Working Directory ---
args <- commandArgs(trailingOnly = TRUE)
WD <- if (length(args) == 0) getwd() else args[1]
setwd(WD)
setwd("D:/GitHub/UrbanDrosophilaEcology")

# --- Data Import and Cleaning ---
DATA <- read.csv("data/Samples_inca_spartacus_vienna_clean_final.csv", header = TRUE)
DATA.Vienna <- na.omit(DATA)
DATA.Vienna <- DATA.Vienna[DATA.Vienna$sampleId != "VCFC_289", ] # Remove Bellaria

# --- Species Data Preparation ---
DATA.spec <- as.data.frame(DATA.Vienna[, 7:19])
rownames(DATA.spec) <- DATA.Vienna$sampleId
DATA.spec[is.na(DATA.spec)] <- 0
DATA.spec.hell <- decostand(DATA.spec, method = "hellinger")
## only keep species present in at least 5 samples
DATA.spec <- DATA.spec[, colSums(DATA.spec) > 5]   


# --- Species Association Analysis ---
library(Hmisc)   # for rcorr
library(tidyverse)
library(igraph)

# Create output folder for species association results
assoc_dir <- "results/species_associations"
dir.create(assoc_dir, showWarnings = FALSE, recursive = TRUE)

# Pairwise Spearman correlations on raw abundance data
assoc_res <- rcorr(as.matrix(DATA.spec), type = "spearman")

# Extract correlations and p-values
corr_mat <- assoc_res$r
pval_mat <- assoc_res$P

# FDR correction
pval_adj_mat <- matrix(
  p.adjust(pval_mat, method = "fdr"), 
  nrow = nrow(pval_mat), 
  dimnames = dimnames(pval_mat)
)

# Save full matrices
write.csv(corr_mat, file.path(assoc_dir, "species_correlations.csv"), quote = FALSE)
write.csv(pval_mat, file.path(assoc_dir, "species_pvalues.csv"), quote = FALSE)
write.csv(pval_adj_mat, file.path(assoc_dir, "species_pvalues_FDR.csv"), quote = FALSE)

# Create a tidy table of significant associations
sig_assoc <- which(pval_adj_mat < 0.01 & lower.tri(pval_adj_mat), arr.ind = TRUE) %>%
  as.data.frame() %>%
  mutate(
    Species1 = rownames(corr_mat)[row],
    Species2 = colnames(corr_mat)[col],
    Correlation = corr_mat[cbind(row, col)],
    p_value = pval_mat[cbind(row, col)],
    p_adj_FDR = pval_adj_mat[cbind(row, col)]
  ) %>%
  select(Species1, Species2, Correlation, p_value, p_adj_FDR) %>%
  arrange(p_adj_FDR)

write.csv(sig_assoc, file.path(assoc_dir, "significant_species_associations.csv"), row.names = FALSE)
## --- Significant Species Association Network ---
if (nrow(sig_assoc) > 0) {
    # Keep only significant associations for the network
    sig_edges <- sig_assoc %>% filter(p_adj_FDR < 0.05)

    # Create igraph object
    g <- graph_from_data_frame(sig_edges)

    # Edge color: green for positive, red for negative
    E(g)$color <- ifelse(E(g)$Correlation > 0, "darkgreen", "red")

    # Edge width proportional to absolute correlation
    E(g)$width <- abs(E(g)$Correlation) * 5

    # Node color based on degree (number of associations)
    degree_vals <- degree(g)
    pal_nodes <- colorRampPalette(c("lightblue", "darkblue"))
    V(g)$color <- pal_nodes(max(degree_vals) + 1)[degree_vals + 1]

    # Save network plot with legends
    pdf(file.path(assoc_dir, "species_association_network_significant.pdf"), width = 9, height = 9)

    plot(
        g,
        vertex.size = 25,
        vertex.label.cex = 0.8,
        main = "Significant Species Associations (FDR < 0.05)"
    )

    # Legend for correlation direction
    legend(
        "topleft",
        legend = c("Positive correlation", "Negative correlation"),
        col = c("darkgreen", "red"),
        lwd = 3,
        bty = "n"
    )

    # Legend for correlation strength (edge width)
    # Legend for correlation strength (edge width) using same scaling as network
    legend_corr_vals <- c(0.1, 0.3, 0.5)
    legend(
        "bottomleft",
        legend = paste0(c("Weak ", "Moderate ", "Strong "), "(", legend_corr_vals, ")"),
        lwd = abs(legend_corr_vals) * 5, # same scaling as E(g)$width
        col = "grey50",
        bty = "n",
        title = "Correlation strength"
    )

    # Legend for node fill colors (degree)
    par(xpd = TRUE) # allow drawing outside plot
    legend_colors <- pal_nodes(5)
    legend_labels <- seq(min(degree_vals), max(degree_vals), length.out = 5)
    legend(
        "topright",
        legend = round(legend_labels, 0),
        fill = legend_colors,
        border = "black",
        title = "Node degree",
        bty = "n"
    )
    par(xpd = FALSE)

    dev.off()
}
