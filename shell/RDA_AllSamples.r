# RDA Analysis for Urban Drosophila Ecology
# Author: Martin Kapun

# Description: This script performs Redundancy Analysis (RDA) on Drosophila sample data,
# including data cleaning, multicollinearity checks, variable selection, scaling, and visualization.

# --- Load Required Libraries ---
library(tidyverse) # Data manipulation and visualization
library(ggmap) # Mapping tools
library(gstat) # Geostatistical modelling
library(osmdata) # For getbb() - OpenStreetMap data
library(FactoMineR) # Multivariate data analysis
library(factoextra) # Visualization for multivariate data
library(vegan) # Ecological analysis (RDA, etc.)
library(ggnewscale) # Multiple fill/color scales in ggplot2
library(ggpubr) # Publication-ready plots
library(corrplot) # Correlation matrix visualization

# --- Set Working Directory ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    WD <- getwd()
} else {
    WD <- args[1]
}
setwd(WD)

# --- Data Import and Preparation ---
# Read cleaned sample data
DATA <- read.csv("data/Samples_inca_spartacus_vienna_clean.csv", header = TRUE)

# Create results directory if it doesn't exist
dir.create("results/RDA_all", showWarnings = FALSE)

# Drop "Viennese" factors (keep columns 1:37)
DATA.all <- DATA[, 1:37]
DATA.all <- na.omit(DATA.all) # Remove rows with missing values

# Select environmental variables (columns 20 to end)
DATA.env.all <- DATA.all %>%
    dplyr::select(20:ncol(DATA.all))

# Remove columns with only zeros
DATA.env.all <- DATA.env.all[, colSums(DATA.env.all) != 0]

# --- Multicollinearity Check ---
cor_matrix <- cor(DATA.env.all, use = "complete.obs")

# Visualize the full correlation matrix
pdf("results/RDA_all/Corrplot_full.pdf", width = 20, height = 20)
corrplot::corrplot(cor_matrix, method = "circle")
dev.off()

# --- Remove Redundant Variables ---
DATA.env.all <- DATA.env.all %>%
    dplyr::select(
        -INCAL_WindSpeedNorth_daily,
        -INCAL_GlobalRadiation_daily,
        -INCAL_Temperature2m_daily,
        -Monthly_TM,
        -INCAL_RainfallRate_daily
    )

# Re-calculate and visualize reduced correlation matrix
cor_matrix <- cor(DATA.env.all, use = "complete.obs")
pdf("results/RDA_all/Corrplot_reduced.pdf", width = 20, height = 20)
corrplot::corrplot(cor_matrix, method = "circle")
dev.off()

# --- Data Scaling ---
DATA.env_scaled.all <- scale(DATA.env.all)

# --- Extract Coordinates, Dates, and Species Counts ---
DATA.coord.all <- DATA.all[, 5:6] # Latitude, Longitude
DATA.date.all <- DATA.all[, 4] # Date

# Species counts (columns 7:19), replace NAs with 0
DATA.spec.all <- DATA.all[, 7:19]
DATA.spec.all[is.na(DATA.spec.all)] <- 0

# Hellinger transformation for species data
DATA.spec.all.hell <- decostand(DATA.spec.all, method = "hellinger")

# Combine environmental, coordinate, and date data for RDA
DATA.factors.all <- cbind(DATA.env.all, DATA.coord.all, DATA.date.all)
DATA.factors.all[is.na(DATA.factors.all)] <- 0

# --- Redundancy Analysis (RDA) ---
# Full model: all variables except Latitude/Longitude, conditioned on Latitude/Longitude
rda_result.all <- rda(
    DATA.spec.all.hell ~ . - Latitude - Longitude + Condition(Latitude + Longitude),
    data = as.data.frame(DATA.factors.all)
)

# Null model: only conditioning on Latitude/Longitude
rda_result.0.all <- rda(
    DATA.spec.all.hell ~ 1 + Condition(Latitude + Longitude),
    data = as.data.frame(DATA.factors.all)
)

# Forward selection of variables using ordiR2step
rda_result.osR2.all <- ordiR2step(
    rda_result.0.all,
    scope = formula(rda_result.all),
    direction = "forward",
    permutations = 99999
)

# --- Output RDA Results ---
sink("results/RDA_all/stats.txt")
summary(rda_result.osR2.all)
anova.cca(rda_result.osR2.all, step = 1000)
anova.cca(rda_result.osR2.all, step = 10000, by = "term")
RsquareAdj(rda_result.osR2.all)
sink()

# --- Ordination Plot ---
pdf("results/RDA_all/Ordiplot.pdf", width = 8, height = 6)
ordiplot(rda_result.osR2.all, scaling = 1, type = "text")
dev.off()
