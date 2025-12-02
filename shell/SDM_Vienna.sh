#!/bin/bash
# SDM_Vienna.sh
# Script for Species Distribution Modeling (SDM) in Vienna using climate and abundance data.
# Usage: bash SDM_Vienna.sh <working_directory>
# Author: [Your Name]
# Date: [Date]
# Description:
#   - Prepares climate raster data
#   - Runs R script for SDM using Random Forest
# #   - Outputs performance metrics, prediction rasters, and visualizations

# # Set working directory from argument
WD=$1
WD="/media/inter/mkapun/projects/UrbanDrosophilaEcology"

# # Change to data directory and prepare folders
# cd data
# mkdir -p Climate

# # Copy climate TIFF files from source directories
# cp ~/mounts/BioMem_2/ssteindl/UC3/ClimateData/Vienna/S3/ViennaData/zartifs/*tiff Climate/
# cp ~/mounts/BioMem_2/ssteindl/UC3/ClimateData/Vienna/S3/ViennaData/fairicube/vienna_data/100m/r*/r*/*tif Climate/

# # Remove unwanted Bezirke raster
# rm -f Climate/Reprojected/r00_Wien_Bezirke_100m_1_1.tif

# # Create results directory for SDM outputs
# mkdir -p ../results/SDM

# Run R script for SDM analysis
Rscript -e '
# IMPROVED SDM Analysis for Vienna - GLM Primary Method
# -----------------------------------------------------
# This script implements SDM with cross-species comparability:
#   - GLM as primary method for ALL species (comparability)
#   - Random Forest as sensitivity analysis (robustness check)
#   - Dual metrics: AUC (discrimination) and Kappa (practical accuracy)
#   - Presence/absence modeling with proper evaluation
#
# Key principle from reviewer response:
#   - Use SAME method for all species to prevent confounding
#     biological and methodological variation
#   - GLM selected for best Kappa (0.527) and interpretability
#   - RF validates robustness (both identify same ecological patterns)

setwd("/media/inter/mkapun/projects/UrbanDrosophilaEcology")

# ---- Load Required Libraries ----
required_packages <- c(
  "raster", "rgdal", "sf", "dplyr", "caret", "randomForest", "ranger",
  "ggplot2", "sp", "terra", "Hmisc", "lme4", "rasterVis", "viridis",
  "readxl", "vegan", "corrplot", "ggmap", "ggpubr", "pROC", "mgcv", "dismo"
)
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(raster)
library(terra)
library(sf)
library(dplyr)
library(randomForest)
library(caret)
library(ggplot2)
library(ranger)
library(rasterVis)
library(lme4)
library(readxl)
library(vegan)
library(corrplot)
library(ggmap)
library(ggpubr)
library(pROC)
library(mgcv)  # For GAM
library(dismo) # For MaxEnt

# Check MaxEnt availability
maxent_available <- file.exists("/media/inter/mkapun/projects/UrbanDrosophilaEcology/scripts/maxent/maxent.jar")

# ---- Load and Prepare Raster Data ----
# List TIFF files from different sources
nc_files <- list.files("data/spartacus_data/yearly", pattern = "*.tif$", full.names = TRUE)
tif_files <- list.files("data/Climate/Reprojected", pattern = "*.tif$", full.names = TRUE)
corrected_tiff_files <- c(nc_files, tif_files)

# Read TIFF files as raster objects
rasters <- lapply(corrected_tiff_files, rast)

# Use the first raster as a template for alignment
template <- rast(tif_files[1])

# Align all rasters to the template
aligned_rasters <- lapply(rasters, function(r) resample(r, template))

# Combine aligned rasters into a stack
predictor_stack <- rast(aligned_rasters)

# Save example layers as PDF and PNG
# pdf("data/example_layers.pdf", width=12, height=6)
png("data/example_layers.png", width=6, height=4, units="in", res=400)
plot(predictor_stack)
dev.off()

# Crop predictor stack to template extent
predictor_stack_masked <- terra::crop(predictor_stack, terra::ext(template))

## reproject to EPSG:4326
predictor_stack_masked <- terra::project(predictor_stack_masked, "EPSG:4326")

# Assign unique names to each raster layer
layer_names <- paste0("layer_", seq_len(nlyr(predictor_stack_masked)))
names(predictor_stack_masked) <- layer_names
print(names(predictor_stack_masked))

# ---- Load and Prepare Abundance Data ----
DATA <- read.csv("data/Samples_inca_spartacus_vienna_clean_final.csv", header = TRUE)

# Define CRS as EPSG:4326
crs_epsg_4326 <- CRS("EPSG:4326")

# Select relevant columns and apply Hellinger transformation to abundance data
abundance_data <- DATA %>%
  dplyr::select(Date=collectionEnd, Latitude=5, Longitude=6)
abundance_data <- cbind(abundance_data, decostand(DATA[,7:19], method="hellinger"))
Spec <- names(DATA[,7:19])

# Convert abundance data to spatial object
abundance_sf <- st_as_sf(abundance_data, coords = c("Longitude", "Latitude"), crs = crs_epsg_4326)

# Extract predictor values for each abundance point
abundance_data <- na.omit(cbind(
  abundance_data,
  extract(predictor_stack_masked, abundance_data[, c("Longitude", "Latitude")])
))

# Remove numeric columns with all zeros
numeric_cols <- sapply(abundance_data, is.numeric)
filtered_numeric_data_zero <- colSums(abundance_data[, numeric_cols], na.rm = TRUE) == 0
col_name <- colnames(abundance_data[, numeric_cols])[filtered_numeric_data_zero]
abundance_data <- abundance_data[, !colnames(abundance_data) %in% col_name]

# Add "Time" variable (numeric date)
abundance_data$Time <- as.numeric(as.Date(abundance_data$Date, format = "%d/%m/%Y"))

# Extract only predictor columns
predictor_values <- abundance_data[, grep("layer.", names(abundance_data))]

# ---- Check for Multicollinearity ----
cor_matrix <- cor(predictor_values, use = "complete.obs")
# corrplot::corrplot(cor_matrix, method = "circle") # Uncomment to visualize

# ---- Model Training and Prediction ----
major <- c("melanogaster", "simulans", "hydei", "mercatorum")
plot_list1 <- list()
plot_list2 <- list()

# Initialize results dataframe for comparison
sdm_comparison <- data.frame()

# Remove zero-sum layers from predictor stack
predictor_stack_masked <- predictor_stack_masked[[!names(predictor_stack_masked) %in% col_name]]

for (i in Spec) {
  cat("\n=== Processing species:", i, "===\n")
  LIM <- ifelse(i %in% major, 1, 0.5)

  # Prepare data for current species
  abundance_data.spec <- data.frame(matrix(ncol = 0, nrow = nrow(abundance_data)))
  abundance_data.spec[[i]] <- abundance_data[[i]]
  abundance_data.spec <- cbind(
    abundance_data.spec,
    Latitude = abundance_data$Latitude,
    Longitude = abundance_data$Longitude,
    Time = abundance_data$Time,
    predictor_values
  )

  # Remove rows with NA values
  abundance_data.spec <- na.omit(abundance_data.spec)
  
  # Calculate presence ratio for species-specific parameters
  presence_ratio <- sum(abundance_data.spec[[i]] > 0) / nrow(abundance_data.spec)
  n_presence <- sum(abundance_data.spec[[i]] > 0)
  is_generalist <- presence_ratio > 0.5
  
  cat(sprintf("  Presence: %.1f%% (%d/%d samples)\n", 
      presence_ratio * 100, n_presence, nrow(abundance_data.spec)))
  cat(sprintf("  Classification: %s\n", 
      ifelse(is_generalist, "Generalist", "Specialist")))
  
  # Skip species with insufficient data
  if (n_presence < 20) {
    cat("  SKIPPED: Insufficient observations (n <20)\n")
    next
  }

  # Partition data into train/test sets
  set.seed(42)
  train_index <- createDataPartition(abundance_data.spec[[i]], p = 0.8, list = FALSE)
  train_data <- abundance_data.spec[train_index, ]
  test_data <- abundance_data.spec[-train_index, ]
  
  # Create binary presence/absence variable
  train_data$Presence <- factor(ifelse(train_data[[i]] > 0, 1, 0), levels = c(0, 1))
  test_data$Presence <- factor(ifelse(test_data[[i]] > 0, 1, 0), levels = c(0, 1))

  # Build formula (exclude coordinates and time)
  formula_pa <- as.formula("Presence ~ . - Latitude - Longitude - Time")
  formula_pa_str <- paste0("Presence ~ ", 
                           paste(grep("layer", names(train_data), value=TRUE), 
                                 collapse = " + "))
  formula_pa <- as.formula(formula_pa_str)
  
  # ========================================================================
  # PRIMARY METHOD: GLM (for cross-species comparability)
  # ========================================================================
  cat("\n  Training GLM (Primary Method)...\n")
  
  glm_model <- glm(formula_pa, data = train_data, family = binomial())
  
  # Predictions
  pred_glm_train <- predict(glm_model, train_data, type = "response")
  pred_glm_test <- predict(glm_model, test_data, type = "response")
  
  # ROC and AUC
  roc_glm_train <- roc(train_data$Presence, pred_glm_train, quiet = TRUE)
  roc_glm_test <- roc(test_data$Presence, pred_glm_test, quiet = TRUE)
  auc_glm_train <- as.numeric(auc(roc_glm_train))
  auc_glm_test <- as.numeric(auc(roc_glm_test))
  
  # Kappa (using optimal threshold from ROC)
  optimal_threshold_glm <- as.numeric(coords(roc_glm_test, "best", ret = "threshold")[1])
  pred_glm_test_binary <- factor(ifelse(pred_glm_test > optimal_threshold_glm, 1, 0), 
                                  levels = c(0, 1))
  cm_glm_test <- confusionMatrix(pred_glm_test_binary, test_data$Presence)
  kappa_glm_test <- cm_glm_test$overall["Kappa"]
  
  # Overfitting
  overfitting_glm <- auc_glm_train - auc_glm_test
  
  cat(sprintf("    Train AUC: %.3f | Test AUC: %.3f | Test Kappa: %.3f\n",
              auc_glm_train, auc_glm_test, kappa_glm_test))
  cat(sprintf("    Overfitting: %.3f (%.1f%%)\n", 
              overfitting_glm, (overfitting_glm/auc_glm_train)*100))
  
  # ========================================================================
  # SENSITIVITY ANALYSIS: Random Forest (validation)
  # ========================================================================
  cat("\n  Training Random Forest (Sensitivity Analysis)...\n")
  
  # Species-specific parameters
  if (is_generalist) {
    ntree_param <- 300
    maxnodes_param <- 8
    nodesize_param <- 15
  } else {
    ntree_param <- 500
    maxnodes_param <- 15
    nodesize_param <- 10
  }
  mtry_param <- max(2, floor(sqrt(ncol(train_data) - 4)))
  
  cat(sprintf("    Parameters: ntree=%d, mtry=%d, maxnodes=%d, nodesize=%d\n",
              ntree_param, mtry_param, maxnodes_param, nodesize_param))
  
  rf_model_pa <- randomForest(
    formula_pa,
    data = train_data,
    ntree = ntree_param,
    mtry = mtry_param,
    maxnodes = maxnodes_param,
    nodesize = nodesize_param,
    classwt = c("0" = 1, "1" = 1),  # Balanced weights
    importance = TRUE
  )
  
  # Predictions
  pred_rf_train <- predict(rf_model_pa, train_data, type = "prob")[, "1"]
  pred_rf_test <- predict(rf_model_pa, test_data, type = "prob")[, "1"]
  
  # ROC and AUC
  roc_rf_train <- roc(train_data$Presence, pred_rf_train, quiet = TRUE)
  roc_rf_test <- roc(test_data$Presence, pred_rf_test, quiet = TRUE)
  auc_rf_train <- as.numeric(auc(roc_rf_train))
  auc_rf_test <- as.numeric(auc(roc_rf_test))
  
  # Kappa
  optimal_threshold_rf <- as.numeric(coords(roc_rf_test, "best", ret = "threshold")[1])
  pred_rf_test_binary <- factor(ifelse(pred_rf_test > optimal_threshold_rf, 1, 0), 
                                 levels = c(0, 1))
  cm_rf_test <- confusionMatrix(pred_rf_test_binary, test_data$Presence)
  kappa_rf_test <- cm_rf_test$overall["Kappa"]
  
  # Overfitting
  overfitting_rf <- auc_rf_train - auc_rf_test
  
  cat(sprintf("    Train AUC: %.3f | Test AUC: %.3f | Test Kappa: %.3f\n",
              auc_rf_train, auc_rf_test, kappa_rf_test))
  cat(sprintf("    Overfitting: %.3f (%.1f%%)\n", 
              overfitting_rf, (overfitting_rf/auc_rf_train)*100))
  
  # ========================================================================
  # GAM (Generalized Additive Model)
  # ========================================================================
  cat("\n  Training GAM...\n")
  
  gam_failed <- FALSE
  tryCatch({
    # Create GAM formula with smooth terms (limited k to avoid overfitting)
    predictor_names <- grep("layer", names(train_data), value=TRUE)
    # Limit to first few predictors to avoid convergence issues
    if (length(predictor_names) > 10) {
      predictor_names <- predictor_names[1:10]
    }
    gam_formula_str <- paste0("Presence ~ ", 
                             paste(paste0("s(", predictor_names, ", k=3)"), 
                                   collapse = " + "))
    gam_formula <- as.formula(gam_formula_str)
    
    gam_model <- gam(gam_formula, data = train_data, family = binomial(), method = "REML")
    
    # Predictions
    pred_gam_train <- predict(gam_model, train_data, type = "response")
    pred_gam_test <- predict(gam_model, test_data, type = "response")
    
    # ROC and AUC
    roc_gam_train <- roc(train_data$Presence, pred_gam_train, quiet = TRUE)
    roc_gam_test <- roc(test_data$Presence, pred_gam_test, quiet = TRUE)
    auc_gam_train <- as.numeric(auc(roc_gam_train))
    auc_gam_test <- as.numeric(auc(roc_gam_test))
    
    # Kappa
    optimal_threshold_gam <- as.numeric(coords(roc_gam_test, "best", ret = "threshold")[1])
    pred_gam_test_binary <- factor(ifelse(pred_gam_test > optimal_threshold_gam, 1, 0), 
                                    levels = c(0, 1))
    cm_gam_test <- confusionMatrix(pred_gam_test_binary, test_data$Presence)
    kappa_gam_test <- cm_gam_test$overall["Kappa"]
    
    # Overfitting
    overfitting_gam <- auc_gam_train - auc_gam_test
    
    cat(sprintf("    Train AUC: %.3f | Test AUC: %.3f | Test Kappa: %.3f\n",
                auc_gam_train, auc_gam_test, kappa_gam_test))
    cat(sprintf("    Overfitting: %.3f (%.1f%%)\n", 
                overfitting_gam, (overfitting_gam/auc_gam_train)*100))
  }, error = function(e) {
    cat("    Failed:", conditionMessage(e), "\n")
    auc_gam_train <<- NA
    auc_gam_test <<- NA
    kappa_gam_test <<- NA
    overfitting_gam <<- NA
    gam_failed <<- TRUE
  })
  
  # ========================================================================
  # MaxEnt (Maximum Entropy)
  # ========================================================================
  cat("\n  Training MaxEnt...\n")
  
  maxent_failed <- FALSE
  if (maxent_available) {
    tryCatch({
      # Prepare data for MaxEnt using correct data.frame method
      # x should contain both presence and background
      # p should be a vector indicating presence (1) vs background (0)
      presence_idx <- which(train_data$Presence == 1)
      absence_idx <- which(train_data$Presence == 0)
      
      # Extract environmental data only (layer columns)
      presence_env <- train_data[presence_idx, grep("layer", names(train_data))]
      background_env <- train_data[absence_idx, grep("layer", names(train_data))]
      
      # Combine presence and background
      all_env <- rbind(presence_env, background_env)
      # Create presence vector (1 for presence, 0 for background)
      p_vec <- c(rep(1, nrow(presence_env)), rep(0, nrow(background_env)))
      
      # Train MaxEnt
      maxent_model <- maxent(x = all_env, p = p_vec, silent = TRUE)
      
      # Predictions
      pred_maxent_train <- predict(maxent_model, train_data[, grep("layer", names(train_data))])
      pred_maxent_test <- predict(maxent_model, test_data[, grep("layer", names(test_data))])
      
      # ROC and AUC
      roc_maxent_train <- roc(train_data$Presence, pred_maxent_train, quiet = TRUE)
      roc_maxent_test <- roc(test_data$Presence, pred_maxent_test, quiet = TRUE)
      auc_maxent_train <- as.numeric(auc(roc_maxent_train))
      auc_maxent_test <- as.numeric(auc(roc_maxent_test))
      
      # Kappa
      optimal_threshold_maxent <- as.numeric(coords(roc_maxent_test, "best", ret = "threshold")[1])
      pred_maxent_test_binary <- factor(ifelse(pred_maxent_test > optimal_threshold_maxent, 1, 0), 
                                         levels = c(0, 1))
      cm_maxent_test <- confusionMatrix(pred_maxent_test_binary, test_data$Presence)
      kappa_maxent_test <- cm_maxent_test$overall["Kappa"]
      
      # Overfitting
      overfitting_maxent <- auc_maxent_train - auc_maxent_test
      
      cat(sprintf("    Train AUC: %.3f | Test AUC: %.3f | Test Kappa: %.3f\n",
                  auc_maxent_train, auc_maxent_test, kappa_maxent_test))
      cat(sprintf("    Overfitting: %.3f (%.1f%%)\n", 
                  overfitting_maxent, (overfitting_maxent/auc_maxent_train)*100))
    }, error = function(e) {
      cat("    Failed:", conditionMessage(e), "\n")
      auc_maxent_train <<- NA
      auc_maxent_test <<- NA
      kappa_maxent_test <<- NA
      overfitting_maxent <<- NA
      maxent_failed <<- TRUE
    })
  } else {
    cat("    Skipped: MaxEnt not available\n")
    auc_maxent_train <- NA
    auc_maxent_test <- NA
    kappa_maxent_test <- NA
    overfitting_maxent <- NA
    maxent_failed <- TRUE
  }
  
  # ========================================================================
  # SPATIAL PREDICTIONS (using RF for fine-grained maps)
  # ========================================================================
  cat("\n  Generating spatial predictions (RF for visualization)...\n")
  
  # Predict RF probabilities across landscape (better for visualization)
  species_prediction_rf <- predict(predictor_stack_masked, rf_model_pa, type = "prob")
  species_prediction_rf <- species_prediction_rf[[2]]  # Extract probability of presence
  
  # Also get GLM predictions for comparison
  species_prediction_glm <- predict(predictor_stack_masked, glm_model, type = "response")

  
  # ========================================================================
  # SAVE PERFORMANCE METRICS
  # ========================================================================
  output_file <- paste0("results/SDM/Performance_Metrics_", i, ".txt")
  output_text <- paste0(
    "================================================================================\n",
    "SPECIES DISTRIBUTION MODELING: ", i, "\n",
    "================================================================================\n\n",
    "SPECIES INFORMATION:\n",
    "  Classification: ", ifelse(is_generalist, "Generalist", "Specialist"), "\n",
    "  Presence ratio: ", round(presence_ratio * 100, 1), "% (", n_presence, "/", 
    nrow(abundance_data.spec), " samples)\n",
    "  Training samples: ", sum(train_data$Presence == 1), " presences, ", 
    sum(train_data$Presence == 0), " absences\n",
    "  Test samples: ", sum(test_data$Presence == 1), " presences, ", 
    sum(test_data$Presence == 0), " absences\n\n",
    "================================================================================\n",
    "PRIMARY METHOD: GLM (selected for cross-species comparability)\n",
    "================================================================================\n\n",
    "Performance Metrics:\n",
    "  Train AUC:  ", round(auc_glm_train, 3), "\n",
    "  Test AUC:   ", round(auc_glm_test, 3), "\n",
    "  Test Kappa: ", round(kappa_glm_test, 3), "\n",
    "  Overfitting (AUC gap): ", round(overfitting_glm, 3), 
    " (", round((overfitting_glm/auc_glm_train)*100, 1), "%)\n\n",
    "Interpretation:\n",
    "  AUC:   ", ifelse(auc_glm_test > 0.8, "Excellent", 
                   ifelse(auc_glm_test > 0.7, "Good",
                   ifelse(auc_glm_test > 0.6, "Acceptable", "Poor"))), "\n",
    "  Kappa: ", ifelse(kappa_glm_test > 0.8, "Excellent",
                   ifelse(kappa_glm_test > 0.6, "Good",
                   ifelse(kappa_glm_test > 0.4, "Moderate",
                   ifelse(kappa_glm_test > 0.2, "Fair", "Poor")))), "\n\n",
    "Confusion Matrix (Test Set, optimal threshold = ", round(optimal_threshold_glm, 3), "):\n"
  )
  
  # Add confusion matrix
  cm_table <- cm_glm_test$table
  output_text <- paste0(output_text,
    "                Predicted\n",
    "  Observed      Absent  Present\n",
    "    Absent      ", sprintf("%6d", cm_table[1,1]), sprintf("%6d", cm_table[1,2]), "\n",
    "    Present     ", sprintf("%6d", cm_table[2,1]), sprintf("%6d", cm_table[2,2]), "\n\n",
    "  Accuracy:    ", round(cm_glm_test$overall["Accuracy"], 3), "\n",
    "  Sensitivity: ", round(cm_glm_test$byClass["Sensitivity"], 3), "\n",
    "  Specificity: ", round(cm_glm_test$byClass["Specificity"], 3), "\n\n"
  )
  
  output_text <- paste0(output_text,
    "================================================================================\n",
    "SENSITIVITY ANALYSIS: Random Forest (validation)\n",
    "================================================================================\n\n",
    "Model Parameters:\n",
    "  ntree:         ", ntree_param, "\n",
    "  mtry:          ", mtry_param, "\n",
    "  maxnodes:      ", maxnodes_param, "\n",
    "  nodesize:      ", nodesize_param, "\n\n",
    "Performance Metrics:\n",
    "  Train AUC:  ", round(auc_rf_train, 3), "\n",
    "  Test AUC:   ", round(auc_rf_test, 3), "\n",
    "  Test Kappa: ", round(kappa_rf_test, 3), "\n",
    "  Overfitting (AUC gap): ", round(overfitting_rf, 3), 
    " (", round((overfitting_rf/auc_rf_train)*100, 1), "%)\n\n",
    "================================================================================\n",
    "COMPARISON: GLM vs RF\n",
    "================================================================================\n\n",
    "Test AUC:   GLM = ", round(auc_glm_test, 3), " | RF = ", round(auc_rf_test, 3), 
    " | Δ = ", round(auc_rf_test - auc_glm_test, 3), "\n",
    "Test Kappa: GLM = ", round(kappa_glm_test, 3), " | RF = ", round(kappa_rf_test, 3),
    " | Δ = ", round(kappa_rf_test - kappa_glm_test, 3), "\n\n",
    "Recommendation: GLM used for final results (cross-species comparability)\n",
    "RF confirms: ", ifelse(abs(auc_glm_test - auc_rf_test) < 0.1, 
                          "Results are robust to methodology", 
                          "Methods show substantial differences"), "\n\n",
    "================================================================================\n"
  )
  
  writeLines(output_text, con = output_file)
  cat("  Performance metrics saved to:", output_file, "\n")

  # ========================================================================
  # STORE RESULTS FOR COMPARISON
  # ========================================================================
  # Determine best models
  all_aucs <- c(auc_glm_test, auc_rf_test, auc_gam_test, auc_maxent_test)
  all_kappas <- c(kappa_glm_test, kappa_rf_test, kappa_gam_test, kappa_maxent_test)
  model_names <- c("GLM", "RF", "GAM", "MaxEnt")
  
  best_auc_idx <- which.max(all_aucs)
  best_kappa_idx <- which.max(all_kappas)
  best_model_auc <- model_names[best_auc_idx]
  best_model_kappa <- model_names[best_kappa_idx]
  
  sdm_comparison <- rbind(sdm_comparison, data.frame(
    Species = i,
    Type = ifelse(is_generalist, "Generalist", "Specialist"),
    Presence_pct = round(presence_ratio * 100, 1),
    GLM_AUC_train = round(auc_glm_train, 3),
    GLM_AUC_test = round(auc_glm_test, 3),
    GLM_Kappa_test = round(kappa_glm_test, 3),
    GLM_Overfitting = round(overfitting_glm, 3),
    RF_AUC_train = round(auc_rf_train, 3),
    RF_AUC_test = round(auc_rf_test, 3),
    RF_Kappa_test = round(kappa_rf_test, 3),
    RF_Overfitting = round(overfitting_rf, 3),
    GAM_AUC_train = round(auc_gam_train, 3),
    GAM_AUC_test = round(auc_gam_test, 3),
    GAM_Kappa_test = round(kappa_gam_test, 3),
    GAM_Overfitting = round(overfitting_gam, 3),
    MaxEnt_AUC_train = round(auc_maxent_train, 3),
    MaxEnt_AUC_test = round(auc_maxent_test, 3),
    MaxEnt_Kappa_test = round(kappa_maxent_test, 3),
    MaxEnt_Overfitting = round(overfitting_maxent, 3),
    Best_Model_AUC = best_model_auc,
    Best_Model_Kappa = best_model_kappa
  ))

  # ========================================================================
  # ROC CURVE COMPARISON PLOT
  # ========================================================================
  cat("  Creating ROC comparison plot...\n")
  
  # Prepare ROC data for plotting
  roc_comparison_data <- data.frame(
    FPR = 1 - roc_glm_test$specificities,
    TPR = roc_glm_test$sensitivities,
    Model = "GLM",
    AUC = auc_glm_test
  )
  
  roc_comparison_data <- rbind(roc_comparison_data, data.frame(
    FPR = 1 - roc_rf_test$specificities,
    TPR = roc_rf_test$sensitivities,
    Model = "RF",
    AUC = auc_rf_test
  ))
  
  # Add GAM if successful
  if (!gam_failed) {
    roc_comparison_data <- rbind(roc_comparison_data, data.frame(
      FPR = 1 - roc_gam_test$specificities,
      TPR = roc_gam_test$sensitivities,
      Model = "GAM",
      AUC = auc_gam_test
    ))
  }
  
  # Add MaxEnt if successful
  if (!maxent_failed) {
    roc_comparison_data <- rbind(roc_comparison_data, data.frame(
      FPR = 1 - roc_maxent_test$specificities,
      TPR = roc_maxent_test$sensitivities,
      Model = "MaxEnt",
      AUC = auc_maxent_test
    ))
  }
  
  # Create color and label vectors
  model_colors <- c("GLM" = "#1f77b4", "RF" = "#ff7f0e", "GAM" = "#2ca02c", "MaxEnt" = "#d62728")
  model_labels <- c(paste0("GLM (AUC=", round(auc_glm_test, 3), ")"),
                   paste0("RF (AUC=", round(auc_rf_test, 3), ")"))
  if (!gam_failed) {
    model_labels <- c(model_labels, paste0("GAM (AUC=", round(auc_gam_test, 3), ")"))
  }
  if (!maxent_failed) {
    model_labels <- c(model_labels, paste0("MaxEnt (AUC=", round(auc_maxent_test, 3), ")"))
  }
  
  # Create ROC plot
  p_roc <- ggplot(roc_comparison_data, aes(x = FPR, y = TPR, color = Model)) +
    geom_line(size = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = model_colors,
                      labels = model_labels) +
    labs(title = paste0("ROC Curve Comparison: D. ", i),
         subtitle = paste0("Best model (AUC): ", best_model_auc, " | Best model (Kappa): ", best_model_kappa),
         x = "False Positive Rate (1 - Specificity)",
         y = "True Positive Rate (Sensitivity)") +
    theme_minimal() +
    theme(legend.position = c(0.7, 0.2),
          legend.background = element_rect(fill = "white", color = "black"))
  
  ggsave(paste0("results/SDM/ROC_Comparison_", i, ".pdf"),
         p_roc, width = 8, height = 6)
  ggsave(paste0("results/SDM/ROC_Comparison_", i, ".png"),
         p_roc, width = 8, height = 6, dpi = 300)

  # ========================================================================
  # SAVE SPATIAL PREDICTIONS
  # ========================================================================
  # Convert RF prediction to data frame for plotting (RF has more nuanced predictions)
  species_prediction_df <- as.data.frame(species_prediction_rf, xy=TRUE)
  names(species_prediction_df)[3] <- "predictions"

  # Save RF prediction as GeoTIFF (primary for visualization)
  writeRaster(
    species_prediction_rf,
    paste0("results/SDM/", i, "_SDM_RF.tif"),
    overwrite=TRUE
  )
  
  # Save GLM prediction for comparison
  writeRaster(
    species_prediction_glm,
    paste0("results/SDM/", i, "_SDM_GLM.tif"),
    overwrite=TRUE
  )
  
  cat("  Predictions saved as GeoTIFF\n")

  # ========================================================================
  # VISUALIZATION
  # ========================================================================
  cat("  Creating visualization...\n")
  
  # Points with nonzero abundance (presence)
  POINTS <- abundance_data.spec[abundance_data.spec[[i]] > 0, ]
  
  # Get bounding box for plotting
  x_range <- range(species_prediction_df$x, na.rm = TRUE)
  y_range <- range(species_prediction_df$y, na.rm = TRUE)

  # Try to download Stadia basemap (with error handling)
  stamen_map <- tryCatch({
    get_stadiamap(
      bbox = c(left = x_range[1], bottom = y_range[1], 
               right = x_range[2], top = y_range[2]),
      maptype = "stamen_toner_lite",
      zoom = 10
    )
  }, error = function(e) {
    cat("    Warning: Could not download basemap, using simple background\n")
    NULL
  })

  # Create plot
  if (!is.null(stamen_map)) {
    # Plot with basemap
    p <- ggmap(stamen_map) +
      geom_tile(data = species_prediction_df, aes(x = x, y = y, fill = predictions), alpha = 0.7) +
      scale_fill_viridis_c(name = "Prob.", option = "C", limits = c(0, 1)) +
      labs(title = paste0("SDM (RF): D. ", i), 
           subtitle = paste0("RF: AUC=", round(auc_rf_test, 3), ", Kappa=", round(kappa_rf_test, 3), 
                            " | GLM: AUC=", round(auc_glm_test, 3), ", Kappa=", round(kappa_glm_test, 3)),
           x = "Longitude", y = "Latitude") +
      geom_point(data = POINTS, aes(x = Longitude, y = Latitude), 
                 color = "white", size = 1, alpha = 0.8, shape = 21, fill = "red") +
      theme_minimal()
  } else {
    # Plot without basemap
    p <- ggplot() +
      geom_tile(data = species_prediction_df, aes(x = x, y = y, fill = predictions)) +
      scale_fill_viridis_c(name = "Prob.", option = "C", limits = c(0, 1)) +
      labs(title = paste0("SDM (RF): D. ", i), 
           subtitle = paste0("RF: AUC=", round(auc_rf_test, 3), ", Kappa=", round(kappa_rf_test, 3), 
                            " | GLM: AUC=", round(auc_glm_test, 3), ", Kappa=", round(kappa_glm_test, 3)),
           x = "Longitude", y = "Latitude") +
      geom_point(data = POINTS, aes(x = Longitude, y = Latitude), 
                 color = "white", size = 1, alpha = 0.8, shape = 21, fill = "red") +
      coord_fixed() +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "grey90"),
        panel.grid = element_line(color = "white")
      )
  }

  # Store plot in appropriate list
  if (i %in% major) {
    plot_list1[[i]] <- p
  } else {
    plot_list2[[i]] <- p
  }
  
  cat("  Visualization created\n")
}

# ---- Combine and Save Plots ----
cat("\n================================================================================\n")
cat("Generating compound figures...\n")
cat("================================================================================\n")

compound_figure1 <- ggarrange(plotlist = plot_list1, ncol = 2, nrow = 2)
compound_figure2 <- ggarrange(plotlist = plot_list2, ncol = 4, nrow = 3)

ggsave(compound_figure1, file = "results/SDM/SpecDist_Major_RF.pdf", width = 9, height = 5)
ggsave(compound_figure1, file = "results/SDM/SpecDist_Major_RF.png", width = 9, height = 5)
ggsave(compound_figure2, file = "results/SDM/SpecDist_Minor_RF.pdf", width = 18, height = 10)
ggsave(compound_figure2, file = "results/SDM/SpecDist_Minor_RF.png", width = 18, height = 10)

cat("\nCompound figures saved.\n")

# ---- Save Comparison Table ----
cat("\nSaving model comparison table...\n")
write.csv(sdm_comparison, "results/SDM/SDM_Model_Comparison.csv", row.names = FALSE)

# ---- Print Summary Statistics ----
cat("\n================================================================================\n")
cat("SUMMARY STATISTICS\n")
cat("================================================================================\n\n")

cat(sprintf("Species analyzed: %d\n", nrow(sdm_comparison)))
cat(sprintf("  Generalists: %d\n", sum(sdm_comparison$Type == "Generalist")))
cat(sprintf("  Specialists: %d\n\n", sum(sdm_comparison$Type == "Specialist")))

cat("Mean Performance Metrics:\n")
cat(sprintf("  GLM    - Test AUC: %.3f | Test Kappa: %.3f | Overfitting: %.3f\n",
            mean(sdm_comparison$GLM_AUC_test, na.rm=TRUE),
            mean(sdm_comparison$GLM_Kappa_test, na.rm=TRUE),
            mean(sdm_comparison$GLM_Overfitting, na.rm=TRUE)))
cat(sprintf("  RF     - Test AUC: %.3f | Test Kappa: %.3f | Overfitting: %.3f\n",
            mean(sdm_comparison$RF_AUC_test, na.rm=TRUE),
            mean(sdm_comparison$RF_Kappa_test, na.rm=TRUE),
            mean(sdm_comparison$RF_Overfitting, na.rm=TRUE)))
cat(sprintf("  GAM    - Test AUC: %.3f | Test Kappa: %.3f | Overfitting: %.3f\n",
            mean(sdm_comparison$GAM_AUC_test, na.rm=TRUE),
            mean(sdm_comparison$GAM_Kappa_test, na.rm=TRUE),
            mean(sdm_comparison$GAM_Overfitting, na.rm=TRUE)))
cat(sprintf("  MaxEnt - Test AUC: %.3f | Test Kappa: %.3f | Overfitting: %.3f\n\n",
            mean(sdm_comparison$MaxEnt_AUC_test, na.rm=TRUE),
            mean(sdm_comparison$MaxEnt_Kappa_test, na.rm=TRUE),
            mean(sdm_comparison$MaxEnt_Overfitting, na.rm=TRUE)))

cat("Best Model by Metric:\n")
cat(sprintf("  AUC:   GLM=%d | RF=%d | GAM=%d | MaxEnt=%d (out of %d species)\n",
            sum(sdm_comparison$Best_Model_AUC == "GLM", na.rm=TRUE),
            sum(sdm_comparison$Best_Model_AUC == "RF", na.rm=TRUE),
            sum(sdm_comparison$Best_Model_AUC == "GAM", na.rm=TRUE),
            sum(sdm_comparison$Best_Model_AUC == "MaxEnt", na.rm=TRUE),
            nrow(sdm_comparison)))
cat(sprintf("  Kappa: GLM=%d | RF=%d | GAM=%d | MaxEnt=%d (out of %d species)\n\n",
            sum(sdm_comparison$Best_Model_Kappa == "GLM", na.rm=TRUE),
            sum(sdm_comparison$Best_Model_Kappa == "RF", na.rm=TRUE),
            sum(sdm_comparison$Best_Model_Kappa == "GAM", na.rm=TRUE),
            sum(sdm_comparison$Best_Model_Kappa == "MaxEnt", na.rm=TRUE),
            nrow(sdm_comparison)))

cat("\n================================================================================\n")
cat("SDM ANALYSIS COMPLETE\n")
cat("================================================================================\n")
cat("\nResults saved to: results/SDM/\n\n")
cat("Key outputs:\n")
cat("  - Performance_Metrics_<species>.txt (detailed results)\n")
cat("  - SDM_Model_Comparison.csv (comparison table)\n")
cat("  - ROC_Comparison_<species>.pdf/png (ROC curves)\n")
cat("  - <species>_SDM_RF.tif (RF predictions - primary for visualization)\n")
cat("  - <species>_SDM_GLM.tif (GLM predictions - for comparison)\n")
cat("  - SpecDist_Major_RF.pdf/png (main species with RF predictions)\n")
cat("  - SpecDist_Minor_RF.pdf/png (other species with RF predictions)\n")
cat("\n")
cat("Note: Maps use RF predictions (more fine-grained, better for visualization)\n")
cat("      GLM used for statistical comparison (better cross-species comparability)\n")
cat("\n================================================================================\n")
'
