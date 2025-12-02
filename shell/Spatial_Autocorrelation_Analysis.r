# Spatial_Autocorrelation_Analysis.r
# Spatial autocorrelation and linear mixed models
# Addressing: "Test for spatial autocorrelation (Moran's I)" and "Use linear mixed models with spatial structure"
# Author: Martin Kapun
# Date: November 2025

cat("\n")
cat("===============================================================\n")
cat("  SPATIAL AUTOCORRELATION ANALYSIS\n")
cat("===============================================================\n\n")

# =============================================================================
# PACKAGE INSTALLATION AND LOADING
# =============================================================================

cat("Checking and installing required packages...\n")

required_packages <- c("tidyverse", "vegan", "spdep", "ncf", "ape", "nlme", "lubridate")

# Function to install missing packages
install_if_missing <- function(packages) {
    new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
    if (length(new_packages) > 0) {
        cat("Installing missing packages:", paste(new_packages, collapse = ", "), "\n")
        install.packages(new_packages, repos = "https://cran.r-project.org", dependencies = TRUE)
    } else {
        cat("All required packages are already installed.\n")
    }
}

install_if_missing(required_packages)

# Load libraries
library(tidyverse)
library(vegan)
library(spdep)
library(ncf)
library(ape)
library(nlme)
library(lubridate)

cat("All libraries loaded successfully!\n\n")

# =============================================================================
# DATA LOADING AND PREPARATION
# =============================================================================

cat("Loading and preparing data...\n")

# Set working directory
WD <- "/media/inter/mkapun/projects/UrbanDrosophilaEcology"
setwd(WD)

# Create output directory
dir.create("results/Spatial_Autocorrelation", recursive = TRUE, showWarnings = FALSE)

# Load data
DATA <- read.csv("data/Samples_inca_spartacus_vienna_clean_final.csv", header = TRUE)
DATA.Vienna <- na.omit(DATA)
DATA.Vienna <- DATA.Vienna[DATA.Vienna$sampleId != "VCFC_289", ] # Remove Bellaria
DATA.Vienna <- DATA.Vienna[rowSums(DATA.Vienna[, 7:19]) > 0, ] # Remove empty samples

cat("Total samples after filtering:", nrow(DATA.Vienna), "\n\n")

# Parse date information
cat("Parsing temporal information...\n")
DATA.Vienna$collectionEnd <- dmy(DATA.Vienna$collectionEnd)
DATA.Vienna$Month <- month(DATA.Vienna$collectionEnd)
DATA.Vienna$Week <- week(DATA.Vienna$collectionEnd)
DATA.Vienna$MonthName <- factor(month.name[DATA.Vienna$Month], levels = month.name)

# Prepare species composition data
cat("Preparing species composition data...\n")
DATA.spec <- as.data.frame(DATA.Vienna[, 7:19])
rownames(DATA.spec) <- DATA.Vienna$sampleId
DATA.spec[is.na(DATA.spec)] <- 0

# Prepare environmental data
cat("Preparing environmental data...\n")
DATA.env.Vienna <- DATA.Vienna %>%
    select(20:ncol(.)) %>%
    select(where(~ is.numeric(.) && sum(., na.rm = TRUE) != 0 || !is.numeric(.)))

DATA.env.Vienna.numeric <- DATA.env.Vienna %>%
    select(where(is.numeric))
DATA.env_scaled.Vienna <- scale(DATA.env.Vienna.numeric)

# Calculate diversity indices
cat("Calculating diversity indices...\n")
shannon_div <- diversity(DATA.spec, index = "shannon")
simpson_div <- diversity(DATA.spec, index = "simpson")
invsimpson_div <- diversity(DATA.spec, index = "invsimpson")
richness <- specnumber(DATA.spec)
evenness <- shannon_div / log(richness)

cat("Data preparation complete.\n\n")

cat("\n=== SPATIAL AUTOCORRELATION ANALYSIS ===\n")
cat("Testing for spatial autocorrelation using Moran's I...\n\n")

# =============================================================================
# 1. MORAN'S I TESTS
# =============================================================================

sink("results/Spatial_Autocorrelation/spatial_autocorrelation.txt")

# Create spatial coordinates
coords <- cbind(DATA.Vienna$Longitude, DATA.Vienna$Latitude)

# Create spatial weights matrix (k-nearest neighbors, k=8)
cat("Creating spatial weights matrix (k=8 nearest neighbors)...\n")
nb <- knn2nb(knearneigh(coords, k = 8))
listw <- nb2listw(nb, style = "W")

# Test Moran's I for diversity indices
cat("\n--- Moran's I: Shannon Diversity ---\n")
moran_shannon <- moran.test(shannon_div, listw, randomisation = TRUE)
print(moran_shannon)

cat("\n--- Moran's I: Species Richness ---\n")
moran_richness <- moran.test(richness, listw, randomisation = TRUE)
print(moran_richness)

cat("\n--- Moran's I: Simpson Diversity ---\n")
moran_simpson <- moran.test(simpson_div, listw, randomisation = TRUE)
print(moran_simpson)

cat("\n--- Moran's I: Inverse Simpson Diversity ---\n")
moran_invsimpson <- moran.test(invsimpson_div, listw, randomisation = TRUE)
print(moran_invsimpson)

cat("\n--- Moran's I: Eveness ---\n")
# Handle NA values for Eveness (occurs when richness = 1)
# Replace NA with mean for Moran's I test or skip if too many NAs
if (sum(is.na(evenness)) / length(evenness) < 0.3) {
    evenness_filled <- evenness
    evenness_filled[is.na(evenness_filled)] <- mean(evenness, na.rm = TRUE)
    moran_eveness <- moran.test(evenness_filled, listw, randomisation = TRUE)
    print(moran_eveness)
    cat("Note:", sum(is.na(evenness)), "NA values replaced with mean for testing.\n")
} else {
    cat("  Too many NA values (", sum(is.na(evenness)), "/", length(evenness), "). Skipping Moran's I test.\n")
}

# Test Moran's I for individual species abundances
cat("\n--- Moran's I: Individual Species ---\n")
for (sp in colnames(DATA.spec)) {
    if (sum(DATA.spec[[sp]] > 0) > 10) { # Only test if enough non-zero values
        moran_sp <- moran.test(DATA.spec[[sp]], listw, randomisation = TRUE)
        cat(paste0("\n", sp, ":\n"))
        cat(paste0(
            "  Moran's I = ", round(moran_sp$estimate[1], 3),
            ", p-value = ", round(moran_sp$p.value, 4), "\n"
        ))
    }
}

# =============================================================================
# 2. MORAN'S I CORRELOGRAM
# =============================================================================

cat("\n--- Moran's I Correlogram ---\n")
# Create distance-based neighbors for correlogram
dists <- as.matrix(dist(coords))
max_dist <- max(dists) * 0.5 # Use half of maximum distance

# Shannon correlogram
cat("Computing distance-based correlogram...\n")
moran_corr_shannon <- correlog(coords[, 1], coords[, 2], shannon_div,
    increment = max_dist / 10, resamp = 999
)
print(summary(moran_corr_shannon))

sink()

# Plot Moran's I correlogram
cat("Saving correlogram plot...\n")
pdf("results/Spatial_Autocorrelation/Morans_I_Correlogram.pdf", width = 10, height = 6)
plot(moran_corr_shannon, main = "Moran's I Correlogram - Shannon Diversity")
dev.off()

# =============================================================================
# 3. SPATIAL LINEAR MIXED MODELS
# =============================================================================

cat("\n=== SPATIAL LINEAR MIXED MODELS ===\n")
cat("Fitting linear mixed-effects models with spatial structure...\n\n")

sink("results/Spatial_Autocorrelation/spatial_lme.txt")

# Prepare data for spatial LME
lme_data <- data.frame(
    Shannon = shannon_div,
    Richness = richness,
    Simpson = simpson_div,
    InvSimpson = invsimpson_div,
    Eveness = evenness,
    Longitude = DATA.Vienna$Longitude,
    Latitude = DATA.Vienna$Latitude,
    Indoors = DATA.Vienna$Indoors,
    Month = DATA.Vienna$Month,
    DATA.env_scaled.Vienna[, 1:min(10, ncol(DATA.env_scaled.Vienna))] # First 10 env vars
)
lme_data <- na.omit(lme_data)

# Add coordinates for spatial correlation
coordinates_df <- lme_data[, c("Longitude", "Latitude")]

# =============================================================================
# 3.1 Shannon Diversity Models
# =============================================================================

cat("\n--- Spatial LME: Shannon Diversity ---\n")

# Basic model without spatial structure
cat("\nFitting basic LME (no spatial structure)...\n")
lme_basic <- lme(Shannon ~ Longitude + Latitude + Indoors + Month,
    random = ~ 1 | Month, data = lme_data, method = "REML"
)

# Model with exponential spatial correlation structure
cat("Attempting spatial LME (exponential correlation)...\n")
lme_spatial <- try(lme(Shannon ~ Longitude + Latitude + Indoors + Month,
    random = ~ 1 | Month,
    correlation = corExp(form = ~ Longitude + Latitude),
    data = lme_data, method = "REML"
), silent = TRUE)

if (!inherits(lme_spatial, "try-error")) {
    cat("\nBasic LME (no spatial structure):\n")
    print(summary(lme_basic))
    cat("\nAIC (basic): ", AIC(lme_basic), "\n")

    cat("\n\nSpatial LME (exponential correlation):\n")
    print(summary(lme_spatial))
    cat("\nAIC (spatial): ", AIC(lme_spatial), "\n")

    # Compare models
    cat("\n\nModel Comparison (LRT):\n")
    print(anova(lme_basic, lme_spatial))
} else {
    cat("\nWarning: Spatial LME model failed to converge. Showing basic LME only.\n")
    cat("\nBasic LME (no spatial structure):\n")
    print(summary(lme_basic))
    cat("\nAIC (basic): ", AIC(lme_basic), "\n")
}

# =============================================================================
# 3.2 Species Richness Models
# =============================================================================

cat("\n\n--- Spatial LME: Species Richness ---\n")

# Basic model
cat("\nFitting basic LME (no spatial structure)...\n")
lme_basic_rich <- lme(Richness ~ Longitude + Latitude + Indoors + Month,
    random = ~ 1 | Month, data = lme_data, method = "REML"
)

# Spatial model
cat("Attempting spatial LME (exponential correlation)...\n")
lme_spatial_rich <- try(lme(Richness ~ Longitude + Latitude + Indoors + Month,
    random = ~ 1 | Month,
    correlation = corExp(form = ~ Longitude + Latitude),
    data = lme_data, method = "REML"
), silent = TRUE)

if (!inherits(lme_spatial_rich, "try-error")) {
    cat("\nBasic LME (no spatial structure):\n")
    print(summary(lme_basic_rich))
    cat("\nAIC (basic): ", AIC(lme_basic_rich), "\n")

    cat("\n\nSpatial LME (exponential correlation):\n")
    print(summary(lme_spatial_rich))
    cat("\nAIC (spatial): ", AIC(lme_spatial_rich), "\n")

    # Compare models
    cat("\n\nModel Comparison (LRT):\n")
    print(anova(lme_basic_rich, lme_spatial_rich))
} else {
    cat("\nWarning: Spatial LME model failed to converge. Showing basic LME only.\n")
    cat("\nBasic LME (no spatial structure):\n")
    print(summary(lme_basic_rich))
    cat("\nAIC (basic): ", AIC(lme_basic_rich), "\n")
}

# =============================================================================
# 3.3 Simpson Diversity Models
# =============================================================================

cat("\n--- Spatial LME: Simpson Diversity ---\n")

# Basic model
cat("\nFitting basic LME (no spatial structure)...\n")
lme_basic_simpson <- lme(Simpson ~ Longitude + Latitude + Indoors + Month,
    random = ~ 1 | Month, data = lme_data, method = "REML"
)

# Spatial model
cat("Attempting spatial LME (exponential correlation)...\n")
lme_spatial_simpson <- try(lme(Simpson ~ Longitude + Latitude + Indoors + Month,
    random = ~ 1 | Month,
    correlation = corExp(form = ~ Longitude + Latitude),
    data = lme_data, method = "REML"
), silent = TRUE)

if (!inherits(lme_spatial_simpson, "try-error")) {
    cat("\nBasic LME (no spatial structure):\n")
    print(summary(lme_basic_simpson))
    cat("\nAIC (basic): ", AIC(lme_basic_simpson), "\n")

    cat("\n\nSpatial LME (exponential correlation):\n")
    print(summary(lme_spatial_simpson))
    cat("\nAIC (spatial): ", AIC(lme_spatial_simpson), "\n")

    cat("\n\nModel Comparison (LRT):\n")
    print(anova(lme_basic_simpson, lme_spatial_simpson))
} else {
    cat("\nWarning: Spatial LME model failed to converge. Showing basic LME only.\n")
    cat("\nBasic LME (no spatial structure):\n")
    print(summary(lme_basic_simpson))
    cat("\nAIC (basic): ", AIC(lme_basic_simpson), "\n")
}

# =============================================================================
# 3.4 Inverse Simpson Diversity Models
# =============================================================================

cat("\n--- Spatial LME: Inverse Simpson Diversity ---\n")

# Basic model
cat("\nFitting basic LME (no spatial structure)...\n")
lme_basic_invsimpson <- lme(InvSimpson ~ Longitude + Latitude + Indoors + Month,
    random = ~ 1 | Month, data = lme_data, method = "REML"
)

# Spatial model
cat("Attempting spatial LME (exponential correlation)...\n")
lme_spatial_invsimpson <- try(lme(InvSimpson ~ Longitude + Latitude + Indoors + Month,
    random = ~ 1 | Month,
    correlation = corExp(form = ~ Longitude + Latitude),
    data = lme_data, method = "REML"
), silent = TRUE)

if (!inherits(lme_spatial_invsimpson, "try-error")) {
    cat("\nBasic LME (no spatial structure):\n")
    print(summary(lme_basic_invsimpson))
    cat("\nAIC (basic): ", AIC(lme_basic_invsimpson), "\n")

    cat("\n\nSpatial LME (exponential correlation):\n")
    print(summary(lme_spatial_invsimpson))
    cat("\nAIC (spatial): ", AIC(lme_spatial_invsimpson), "\n")

    cat("\n\nModel Comparison (LRT):\n")
    print(anova(lme_basic_invsimpson, lme_spatial_invsimpson))
} else {
    cat("\nWarning: Spatial LME model failed to converge. Showing basic LME only.\n")
    cat("\nBasic LME (no spatial structure):\n")
    print(summary(lme_basic_invsimpson))
    cat("\nAIC (basic): ", AIC(lme_basic_invsimpson), "\n")
}

# =============================================================================
# 3.5 Eveness Models
# =============================================================================

cat("\n--- Spatial LME: Eveness ---\n")

# Basic model
cat("\nFitting basic LME (no spatial structure)...\n")
lme_basic_eveness <- lme(Eveness ~ Longitude + Latitude + Indoors + Month,
    random = ~ 1 | Month, data = lme_data, method = "REML"
)

# Spatial model
cat("Attempting spatial LME (exponential correlation)...\n")
lme_spatial_eveness <- try(lme(Eveness ~ Longitude + Latitude + Indoors + Month,
    random = ~ 1 | Month,
    correlation = corExp(form = ~ Longitude + Latitude),
    data = lme_data, method = "REML"
), silent = TRUE)

if (!inherits(lme_spatial_eveness, "try-error")) {
    cat("\nBasic LME (no spatial structure):\n")
    print(summary(lme_basic_eveness))
    cat("\nAIC (basic): ", AIC(lme_basic_eveness), "\n")

    cat("\n\nSpatial LME (exponential correlation):\n")
    print(summary(lme_spatial_eveness))
    cat("\nAIC (spatial): ", AIC(lme_spatial_eveness), "\n")

    cat("\n\nModel Comparison (LRT):\n")
    print(anova(lme_basic_eveness, lme_spatial_eveness))
} else {
    cat("\nWarning: Spatial LME model failed to converge. Showing basic LME only.\n")
    cat("\nBasic LME (no spatial structure):\n")
    print(summary(lme_basic_eveness))
    cat("\nAIC (basic): ", AIC(lme_basic_eveness), "\n")
}

sink()

cat("\n=== SPATIAL AUTOCORRELATION ANALYSIS COMPLETE ===\n")
cat("Results saved to:\n")
cat("  - results/Spatial_Autocorrelation/spatial_autocorrelation.txt\n")
cat("  - results/Spatial_Autocorrelation/spatial_lme.txt\n")
cat("  - results/Spatial_Autocorrelation/Morans_I_Correlogram.pdf\n\n")
