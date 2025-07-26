#!/bin/bash

###############################################################################
# Data Analysis Pipeline for the ecological analysis of the Vienna City Fly Project
# This script orchestrates the analysis of biodiversity data, including:
#   - Downloading and linking Earth Observation data
#   - Generating descriptive statistics and plots
#   - Performing biodiversity and redundancy analyses
#   - Running species distribution models
#
# Usage:
#   bash main.sh
#
# Requirements:
#   - R and required R packages
#   - Bash shell
#   - All referenced scripts must be present and executable
###############################################################################

# Set working directory
WD="/media/inter/mkapun/projects/UrbanDrosophilaEcology"

echo "Starting Data Analysis Pipeline..."
echo "Working directory: $WD"
echo

# Step 1: Obtain Earth Observation data and link to sampling sites
echo "Step 1: Obtaining Earth Observation data..."
bash "${WD}/shell/GetEOdata.sh"
echo "Earth Observation data obtained."
echo

# Note: Sampling location information (Indoors/Outdoors) was manually added to:
#   Samples_inca_spartacus_vienna_clean_final.csv

# Step 2: Produce plots with descriptive statistics on abundance and temporal patterns
echo "Step 2: Generating descriptive statistics and plots..."
Rscript "${WD}/shell/descriptive.r" "${WD}"
echo "Descriptive statistics and plots generated."
echo

# Step 3: Analyze biodiversity for Viennese samples only
echo "Step 3: Biodiversity analysis for Viennese samples..."
Rscript "${WD}/shell/BioDiv_Vienna.r" "${WD}"
echo "Biodiversity analysis completed."
echo

# Step 4: Redundancy Analyses (RDA) for all samples
echo "Step 4: Redundancy Analysis for all samples..."
Rscript "${WD}/shell/RDA_AllSamples.r" "${WD}"
echo "RDA for all samples completed."
echo

# Step 5: Redundancy Analyses for Vienna samples only
echo "Step 5: Redundancy Analysis for Vienna samples..."
Rscript "${WD}/shell/RDA_Vienna_full.r" "${WD}"
echo "RDA for Vienna samples completed."
echo

# Step 6: Redundancy Analyses with all sampling dates collapsed
echo "Step 6: RDA with collapsed sampling dates..."
Rscript "${WD}/shell/RDA_Vienna_full_collapsed.r" "${WD}"
echo "RDA with collapsed dates completed."
echo

# Step 7: Species Distribution Models
echo "Step 7: Running Species Distribution Models..."
bash "${WD}/shell/SDM_Vienna.sh" "${WD}"
echo "Species Distribution Models completed."
echo

echo "Data Analysis Pipeline finished successfully."
