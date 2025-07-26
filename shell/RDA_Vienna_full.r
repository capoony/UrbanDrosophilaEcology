# RDA Analysis for Vienna Drosophila Data
# Author: Martin Kapun
# ----------------------------------------
# This script performs Redundancy Analysis (RDA) on Drosophila ecology data from Vienna.
# It includes data loading, cleaning, multicollinearity checks, variable selection,
# Hellinger transformation, RDA modeling, and result visualization.

# --- Load Required Packages ---
library(tidyverse) # Data manipulation and visualization
library(vegan) # Ecological analysis (RDA, Hellinger, etc.)
library(corrplot) # Correlation matrix visualization

# --- Set Working Directory ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    WD <- getwd()
} else {
    WD <- args[1]
}
setwd(WD)

# --- Create Output Directory ---
dir.create("results/RDA_Vienna_full", showWarnings = FALSE, recursive = TRUE)

# --- Load and Prepare Data ---
DATA <- read.csv("data/Samples_inca_spartacus_vienna_clean_final.csv", header = TRUE)
DATA.Vienna <- na.omit(DATA) # Remove rows with missing values

# --- Select Environmental Variables ---
DATA.env.Vienna <- DATA.Vienna %>%
    dplyr::select(20:ncol(DATA.Vienna)) # Select environmental columns

DATA.env.Vienna$PI <- DATA.Vienna$ParticipantId # Add ParticipantId

# Remove columns with only zeros (excluding PI)
DATA.env.Vienna <- DATA.env.Vienna[, colSums(DATA.env.Vienna[, 1:(ncol(DATA.env.Vienna) - 1)]) != 0 | names(DATA.env.Vienna) == "PI"]

# --- Multicollinearity Check ---
cor_matrix <- cor(DATA.env.Vienna, use = "complete.obs")

# Visualize the correlation matrix (full)
pdf("results/RDA_Vienna_full/Corrplot_full.pdf", width = 20, height = 20)
corrplot::corrplot(cor_matrix, method = "circle")
dev.off()

# --- Remove Redundant Variables (based on correlation or prior knowledge) ---
DATA.env.Vienna <- DATA.env.Vienna %>%
    dplyr::select(
        -INCAL_WindSpeedNorth_daily,
        -INCAL_GlobalRadiation_daily,
        -INCAL_Temperature2m_daily,
        -Monthly_TM,
        -INCAL_RainfallRate_daily,
        -outdoorsports,
        -vineyard,
        -water.1,
        -commercial
    )

# Re-check and visualize reduced correlation matrix
cor_matrix <- cor(DATA.env.Vienna, use = "complete.obs")
pdf("results/RDA_Vienna_full/Corrplot_reduced.pdf", width = 20, height = 20)
corrplot::corrplot(cor_matrix, method = "circle")
dev.off()

# --- Scale Environmental Variables ---
DATA.env_scaled.Vienna <- scale(DATA.env.Vienna[, 1:(ncol(DATA.env.Vienna) - 1)])
DATA.env_scaled.Vienna <- as.data.frame(DATA.env_scaled.Vienna)
DATA.env_scaled.Vienna$PI <- DATA.env.Vienna$PI

# --- Extract Coordinates and Dates ---
DATA.coord.Vienna <- DATA.Vienna[, 5:6] # Latitude, Longitude
DATA.date.Vienna <- DATA.Vienna[, 4] # Date

# --- Prepare Species Data and Hellinger Transform ---
DATA.spec.Vienna <- DATA.Vienna[, 7:19]
DATA.spec.Vienna[is.na(DATA.spec.Vienna)] <- 0
DATA.spec.Vienna.hell <- decostand(DATA.spec.Vienna, method = "hellinger")

# --- Combine Factors for RDA ---
DATA.factors.Vienna <- cbind(DATA.env.Vienna, DATA.coord.Vienna, DATA.date.Vienna)
DATA.factors.Vienna[is.na(DATA.factors.Vienna)] <- 0

# --- Run RDA Analysis ---
# Full model: all variables except Latitude, Longitude, PI (PI as condition)
rda_result.Vienna <- rda(
    DATA.spec.Vienna.hell ~ . - Latitude - Longitude - PI + Condition(PI),
    data = as.data.frame(DATA.factors.Vienna)
)

# Null model: only PI as condition
rda_result.0.Vienna <- rda(
    DATA.spec.Vienna.hell ~ 1 + Condition(PI),
    data = as.data.frame(DATA.factors.Vienna)
)

# Forward selection of variables
rda_result.osR2.Vienna <- ordiR2step(
    rda_result.0.Vienna,
    scope = formula(rda_result.Vienna),
    direction = "forward",
    permutations = 99999
)

# --- Output RDA Results ---
sink("results/RDA_Vienna_full/stats.txt")
summary(rda_result.osR2.Vienna)
anova.cca(rda_result.osR2.Vienna, step = 1000)
anova.cca(rda_result.osR2.Vienna, step = 10000, by = "term")
RsquareAdj(rda_result.osR2.Vienna)
sink()

# --- Visualize RDA Ordination ---
pdf("results/RDA_Vienna_full/Ordiplot.pdf", width = 8, height = 6)
ordiplot(rda_result.osR2.Vienna, scaling = 1, type = "text")
dev.off()
