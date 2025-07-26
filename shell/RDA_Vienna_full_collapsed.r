# RDA Analysis for Vienna Full Collapsed Dataset
# Author: Martin Kapun
#
# This script performs Redundancy Analysis (RDA) on the Vienna dataset.
# It includes data loading, cleaning, variable selection, multicollinearity checks,
# Hellinger transformation, RDA modeling, and result visualization.

# --- Load Required Libraries ---
library(tidyverse)
library(vegan) # For RDA and ecological analysis
library(corrplot) # For correlation matrix visualization

# --- Set Working Directory ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    WD <- getwd()
} else {
    WD <- args[1]
}
setwd(WD)

# --- Create Results Directory ---
dir.create("results/RDA_Vienna_full_collapsed", showWarnings = FALSE, recursive = TRUE)

# --- Load and Prepare Data ---
DATA <- read.csv("data/Samples_inca_spartacus_vienna_clean_final.csv", header = TRUE)

# Remove rows with missing values
DATA.Vienna <- na.omit(DATA)

# --- Environmental Variables: Selection and Aggregation ---
# Select relevant columns and aggregate by ParticipantId and Indoors
DATA.env.Vienna <- DATA.Vienna %>%
    dplyr::select(2, 20, 29:ncol(DATA.Vienna)) %>%
    group_by(ParticipantId, Indoors) %>%
    summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
    dplyr::select(-ParticipantId)

# Remove columns with only zeros
DATA.env.Vienna <- DATA.env.Vienna[, colSums(DATA.env.Vienna) != 0]

# --- Multicollinearity Check: Correlation Matrix ---
cor_matrix <- cor(DATA.env.Vienna, use = "complete.obs")

# Visualize the full correlation matrix
pdf("results/RDA_Vienna_full_collapsed/Corrplot_full.pdf", width = 20, height = 20)
corrplot(cor_matrix, method = "circle")
dev.off()

# --- Remove Redundant Variables (based on correlation or prior knowledge) ---
DATA.env.Vienna <- DATA.env.Vienna %>%
    dplyr::select(
        -Monthly_TM,
        -Monthly_RR,
        -Monthly_SA,
        -Daily_TN,
        -Daily_RR,
        -Daily_SA,
        -outdoorsports,
        -vineyard,
        -water.1,
        -commercial
    )

# Re-calculate and visualize the reduced correlation matrix
cor_matrix <- cor(DATA.env.Vienna, use = "complete.obs")
pdf("results/RDA_Vienna_full_collapsed/Corrplot_reduced.pdf", width = 20, height = 20)
corrplot(cor_matrix, method = "circle")
dev.off()

# --- Scale Environmental Variables ---
DATA.env_scaled.Vienna <- scale(DATA.env.Vienna)

# --- Extract Coordinates and Dates (if needed for further analysis) ---
DATA.coord.Vienna <- DATA.Vienna[, 5:6]
DATA.date.Vienna <- DATA.Vienna[, 4]

# --- Species Data Preparation and Hellinger Transformation ---
# Aggregate species counts by ParticipantId and Indoors
DATA.spec.Vienna <- DATA.Vienna %>%
    dplyr::select(2, 7:20) %>%
    group_by(ParticipantId, Indoors) %>%
    summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
    dplyr::select(-ParticipantId, -Indoors)

# Replace NAs with zeros
DATA.spec.Vienna[is.na(DATA.spec.Vienna)] <- 0

# Apply Hellinger transformation
DATA.spec.Vienna.hell <- decostand(DATA.spec.Vienna, method = "hellinger")

# --- Redundancy Analysis (RDA) ---
# Full model with all environmental variables
rda_result.Vienna <- rda(DATA.spec.Vienna.hell ~ ., data = as.data.frame(DATA.env.Vienna))

# Null model (intercept only)
rda_result.0.Vienna <- rda(DATA.spec.Vienna.hell ~ 1, data = as.data.frame(DATA.env.Vienna))

# Forward selection of variables using ordiR2step
rda_result.osR2.Vienna <- ordiR2step(
    rda_result.0.Vienna,
    scope = formula(rda_result.Vienna),
    direction = "forward",
    permutations = 99999
)

# --- Output RDA Results ---
sink("results/RDA_Vienna_full_collapsed/stats.txt")
summary(rda_result.osR2.Vienna)
anova.cca(rda_result.osR2.Vienna, step = 1000)
anova.cca(rda_result.osR2.Vienna, step = 10000, by = "term")
RsquareAdj(rda_result.osR2.Vienna)
sink()

# --- RDA Ordination Plot ---
pdf("results/RDA_Vienna_full_collapsed/Ordiplot.pdf", width = 8, height = 6)
ordiplot(rda_result.osR2.Vienna, scaling = 1, type = "text")
dev.off()
