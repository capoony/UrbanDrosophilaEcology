# BioDiv_Vienna.r
# Analysis of Drosophila biodiversity in Vienna
# Author: [Your Name]
# Date: [Date]
# Description: This script performs PCA, diversity estimation, NMDS ordination, and mixed-effects modeling on Vienna Drosophila samples.

# --- Load Required Libraries ---
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(vegan) # Biodiversity estimators
library(gridExtra)
library(lme4)
library(car)
library(ggplot2)
library(reshape2)
library(scales)

# --- Set Working Directory ---
args <- commandArgs(trailingOnly = TRUE)
WD <- if (length(args) == 0) getwd() else args[1]
setwd(WD)

# --- Data Import and Cleaning ---
DATA <- read.csv("data/Samples_inca_spartacus_vienna_clean_final.csv", header = TRUE)
DATA.Vienna <- na.omit(DATA)
DATA.Vienna <- DATA.Vienna[DATA.Vienna$sampleId != "VCFC_289", ] # Remove Bellaria
dir.create("results/BioDiv_Vienna", showWarnings = FALSE)
DATA.Vienna <- DATA.Vienna[rowSums(DATA.Vienna[, 7:19]) > 0, ] # Remove empty samples

# --- Environmental Data Preparation ---
DATA.env.Vienna <- DATA.Vienna %>%
    select(20:ncol(.)) %>%
    select(where(~ sum(.) != 0)) # Remove columns with only zeros
DATA.env_scaled.Vienna <- scale(DATA.env.Vienna)

# --- PCA Analysis ---
pca_result <- PCA(DATA.env_scaled.Vienna, scale.unit = FALSE)
scree_plot <- fviz_screeplot(pca_result, addlabels = TRUE, ylim = c(0, 50))
ggsave("results/BioDiv_Vienna/PCA_scree.pdf", scree_plot, width = 10, height = 6)

# Biplots for PCA axes
biplot_1_2 <- fviz_pca_biplot(pca_result, axes = c(1, 2), geom.ind = "blank", repel = TRUE, col.var = "black")
biplot_3_4 <- fviz_pca_biplot(pca_result, axes = c(3, 4), geom.ind = "blank", repel = TRUE, col.var = "black")
biplot <- grid.arrange(biplot_1_2, biplot_3_4, ncol = 2)
ggsave("results/BioDiv_Vienna/PCA_biplot.pdf", biplot, width = 15, height = 7)
ggsave("results/BioDiv_Vienna/PCA_biplot.png", biplot, width = 15, height = 7)

# Save PCA loadings
write.table(pca_result$var$coord, "results/BioDiv_Vienna/PCA_loadings.csv", row.names = TRUE, quote = FALSE, sep = ",")

# --- Species Data Preparation ---
DATA.spec <- as.data.frame(DATA.Vienna[, 7:19])
rownames(DATA.spec) <- DATA.Vienna$sampleId
DATA.spec[is.na(DATA.spec)] <- 0
DATA.spec.hell <- decostand(DATA.spec, method = "hellinger")

# --- Diversity Indices Calculation ---
shannon_div <- diversity(DATA.spec, index = "shannon")
simpson_div <- diversity(DATA.spec, index = "simpson")
invsimpson_div <- diversity(DATA.spec, index = "invsimpson")
richness <- specnumber(DATA.spec)
evenness <- shannon_div / log(richness)

# --- Combine Factors and PCA Scores ---
pca_scores <- scale(pca_result$ind$coord[, 1:4])
DATA.factors <- cbind(DATA.Vienna[, 5:6], Days = DATA.Vienna$Days, pca_scores)
DATA.factors[is.na(DATA.factors)] <- 0

# --- Compile Diversity Data ---
div_data <- data.frame(
    Shannon = shannon_div,
    Simpson = simpson_div,
    InvSimpson = invsimpson_div,
    Richness = richness,
    Eveness = evenness,
    DATA.factors
)
write.table(div_data, "results/BioDiv_Vienna/Diversity.csv", quote = FALSE, sep = ",", row.names = TRUE)

# --- NMDS Ordination ---
nmds <- metaMDS(DATA.spec, distance = "bray", k = 2, trymax = 100)
pdf("results/BioDiv_Vienna/Bray_NMDS.pdf", width = 10, height = 6)
plot(nmds, main = "NMDS - Bray-Curtis")
abline(h = 0, v = 0, lty = 2, col = "gray")
text(nmds, display = "species", labels = colnames(DATA.spec), col = "red", cex = 0.8, pos = 3)

# --- Environmental Fit ---
ef <- envfit(nmds, DATA.factors[, 4:7], permutations = 99999)
sink("results/BioDiv_Vienna/Envfit.txt")
print(ef)
sink()
plot(ef, add = TRUE)
dev.off()

# --- Outlier Removal Function ---
remove_outliers <- function(df, cols) {
    for (col in cols) {
        Q1 <- quantile(df[[col]], 0.25, na.rm = TRUE)
        Q3 <- quantile(df[[col]], 0.75, na.rm = TRUE)
        IQR <- Q3 - Q1
        lower <- Q1 - 1.5 * IQR
        upper <- Q3 + 1.5 * IQR
        df <- df %>% filter(.data[[col]] >= lower & .data[[col]] <= upper)
    }
    df
}

# --- Mixed-Effects Modeling and Plotting ---
plot_fixed_effects <- c("Dim.1", "Dim.2", "Dim.3", "Dim.4")
dir.create("results/BioDiv_Vienna/plots", showWarnings = FALSE)
sink("results/BioDiv_Vienna/stats.txt")

for (resp in names(div_data)[1:5]) {
    cat(resp, "\n____________\n")
    formula <- as.formula(paste0(resp, " ~ Latitude + Longitude + Dim.1 + Dim.2 + Dim.3 + Dim.4 + (1 | Days)"))
    cleaned_data <- div_data %>%
        select(all_of(c(resp, "Latitude", "Longitude", plot_fixed_effects, "Days"))) %>%
        remove_outliers(cols = c(resp, "Latitude", "Longitude", plot_fixed_effects))
    TEST <- try(lmer(formula, data = cleaned_data), silent = TRUE)
    if (!inherits(TEST, "try-error")) {
        anova_results <- Anova(TEST, type = "III")
        print(anova_results)
        fixed_effects_p <- data.frame(
            FixedEffect = rownames(anova_results),
            p_value = anova_results$`Pr(>Chisq)`
        ) %>% filter(!is.na(p_value))
    } else {
        fixed_effects_p <- data.frame(FixedEffect = plot_fixed_effects, p_value = NA)
    }
    plot_data_long <- melt(cleaned_data, id.vars = resp, variable.name = "FixedEffect", value.name = "Value") %>%
        filter(FixedEffect %in% plot_fixed_effects) %>%
        left_join(fixed_effects_p, by = "FixedEffect")
    p <- ggplot(plot_data_long, aes(x = Value, y = .data[[resp]])) +
        geom_point(alpha = 0.3) +
        geom_smooth(method = "lm", formula = y ~ x, color = "blue") +
        facet_wrap(~FixedEffect, scales = "free", nrow = 1) +
        labs(title = resp, x = "Fixed Effect Value", y = resp) +
        theme_bw() +
        geom_text(aes(x = Inf, y = Inf, label = paste0("p = ", signif(p_value, 3))),
            inherit.aes = FALSE, hjust = 1.2, vjust = 1.5, size = 3, color = "red"
        ) +
        scale_y_continuous(labels = label_number(accuracy = 0.1))
    ggsave(paste0("results/BioDiv_Vienna/plots/", resp, "_correlations_Dim1-3.png"), p, width = 6, height = 2)
    ggsave(paste0("results/BioDiv_Vienna/plots/", resp, "_correlations_Dim1-3.pdf"), p, width = 6, height = 2)
}
sink()
