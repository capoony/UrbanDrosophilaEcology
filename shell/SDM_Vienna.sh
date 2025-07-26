#!/bin/bash
# SDM_Vienna.sh
# Script for Species Distribution Modeling (SDM) in Vienna using climate and abundance data.
# Usage: bash SDM_Vienna.sh <working_directory>
# Author: [Your Name]
# Date: [Date]
# Description:
#   - Prepares climate raster data
#   - Runs R script for SDM using Random Forest
#   - Outputs performance metrics, prediction rasters, and visualizations

# Set working directory from argument
WD=$1

# Change to data directory and prepare folders
cd data
mkdir -p Climate

# Copy climate TIFF files from source directories
cp ~/mounts/BioMem_2/ssteindl/UC3/ClimateData/Vienna/S3/ViennaData/zartifs/*tiff Climate/
cp ~/mounts/BioMem_2/ssteindl/UC3/ClimateData/Vienna/S3/ViennaData/fairicube/vienna_data/100m/r*/r*/*tif Climate/

# Remove unwanted Bezirke raster
rm -f Climate/Reprojected/r00_Wien_Bezirke_100m_1_1.tif

# Create results directory for SDM outputs
mkdir -p ../results/SDM

# Run R script for SDM analysis
Rscript -e '
# SDM Analysis for Vienna
# -----------------------
# This script:
#   - Loads and processes climate raster data
#   - Loads abundance data
#   - Trains Random Forest models for each species
#   - Evaluates model performance
#   - Predicts and visualizes species distributions

setwd("'"${WD}"'")

# ---- Load Required Libraries ----
required_packages <- c(
  "raster", "rgdal", "sf", "dplyr", "caret", "randomForest", "ranger",
  "ggplot2", "sp", "terra", "Hmisc", "lme4", "rasterVis", "viridis",
  "readxl", "vegan", "corrplot", "ggmap", "ggpubr"
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
pdf("data/example_layers.pdf", width=12, height=6)
png("data/example_layers.png", width=6, height=4, units="in", res=400)
plot(predictor_stack)
dev.off()

# Crop predictor stack to template extent
predictor_stack_masked <- terra::crop(predictor_stack, terra::ext(template))

# Assign unique names to each raster layer
layer_names <- paste0("layer_", seq_len(nlyr(predictor_stack_masked)))
names(predictor_stack_masked) <- layer_names
print(names(predictor_stack_masked))

# ---- Load and Prepare Abundance Data ----
DATA <- read.csv("data/Samples_inca_spartacus_vienna_clean.csv", header = TRUE)

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
abundance_data$Time <- as.numeric(as.Date(abundance_data$Date, format = "%Y-%m-%d"))

# Extract only predictor columns
predictor_values <- abundance_data[, grep("layer.", names(abundance_data))]

# ---- Check for Multicollinearity ----
cor_matrix <- cor(predictor_values, use = "complete.obs")
# corrplot::corrplot(cor_matrix, method = "circle") # Uncomment to visualize

# ---- Model Training and Prediction ----
major <- c("melanogaster", "simulans", "hydei", "mercatorum")
plot_list1 <- list()
plot_list2 <- list()

# Remove zero-sum layers from predictor stack
predictor_stack_masked <- predictor_stack_masked[[!names(predictor_stack_masked) %in% col_name]]

for (i in Spec) {
  cat("Processing species:", i, "\n")
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

  # Partition data into train/test sets
  set.seed(42)
  train_index <- createDataPartition(abundance_data.spec[[i]], p = 0.8, list = FALSE)
  train_data <- abundance_data.spec[train_index, ]
  test_data <- abundance_data.spec[-train_index, ]

  # Build Random Forest formula
  formula <- as.formula(paste0(i, " ~ . - Latitude - Longitude - Time"))

  # Train Random Forest model
  rf_model <- ranger(
    formula = formula,
    data = train_data,
    mtry = 3,
    num.trees = 500,
    importance = "permutation"
  )

  # Predict across raster stack
  species_prediction <- predict(predictor_stack_masked, rf_model, type = "response")

  # Evaluate model on training data
  train_predictions <- predict(rf_model, train_data)$predictions
  train_observed <- train_data[[i]]
  r_squared_train <- cor(train_observed, train_predictions)^2
  rmse_train <- sqrt(mean((train_observed - train_predictions)^2))

  # Evaluate model on test data
  test_predictions <- predict(rf_model, test_data)$predictions
  test_observed <- test_data[[i]]
  r_squared_test <- cor(test_observed, test_predictions)^2
  rmse_test <- sqrt(mean((test_observed - test_predictions)^2))
  mae_test <- mean(abs(test_observed - test_predictions))

  # 5-fold cross-validation
  cv_control <- trainControl(method = "cv", number = 5)
  cv_model <- train(
    formula,
    data = abundance_data.spec,
    method = "ranger",
    trControl = cv_control,
    tuneGrid = expand.grid(mtry = 3, splitrule = "variance", min.node.size = 5)
  )
  cv_results <- cv_model$results
  cv_r_squared <- cv_results$Rsquared[1]
  cv_rmse <- cv_results$RMSE[1]

  # Write performance metrics to file
  output_file <- paste0("results/SDM/Performance_Metrics_", i, ".txt")
  output_text <- paste(
    "Performance Metrics for Species:", i, "\n",
    "Training R-squared:", round(r_squared_train, 3), "\n",
    "Training RMSE:", round(rmse_train, 3), "\n",
    "Test R-squared:", round(r_squared_test, 3), "\n",
    "Test RMSE:", round(rmse_test, 3), "\n",
    "Test MAE:", round(mae_test, 3), "\n",
    sep = ""
  )
  writeLines(output_text, con = output_file)
  cat("Performance metrics saved to:", output_file, "\n")

  # Convert prediction raster to data frame for plotting
  species_prediction_df <- as.data.frame(species_prediction, xy=TRUE)
  names(species_prediction_df)[3] <- "predictions"

  # Save prediction as GeoTIFF
  writeRaster(
    species_prediction,
    paste0("results/SDM/", i, "_SDM.tif"),
    overwrite=TRUE
  )

  # Get bounding box for plotting
  x_range <- range(species_prediction_df$x)
  y_range <- range(species_prediction_df$y)

  # Download Stadia basemap for the bounding box
  stamen_map <- get_stadiamap(
    bbox = c(left = x_range[1], bottom = y_range[1], right = x_range[2], top = y_range[2]),
    maptype = "stamen_toner_lite",
    zoom = 10
  )

  # Points with nonzero abundance
  POINTS <- abundance_data.spec[abundance_data.spec[[i]] != 0, ]

  # Plot prediction on basemap
  p <- ggmap(stamen_map) +
    geom_tile(data = species_prediction_df, aes(x = x, y = y, fill = predictions), alpha = 0.7) +
    scale_fill_viridis_c(name = "Abundance", option = "C", limits = c(0, LIM)) +
    labs(title = paste("SDM for D.", i), x = "Longitude", y = "Latitude") +
    geom_point(data = POINTS, aes(x = Longitude, y = Latitude), color = "white", size = 0.4, alpha = 0.5) +
    theme_minimal()

  # Store plot in appropriate list
  if (i %in% major) {
    plot_list1[[i]] <- p
  } else {
    plot_list2[[i]] <- p
  }
}

# ---- Combine and Save Plots ----
compound_figure1 <- ggarrange(plotlist = plot_list1, ncol = 2, nrow = 2)
compound_figure2 <- ggarrange(plotlist = plot_list2, ncol = 4, nrow = 3)

ggsave(compound_figure1, file = "results/SDM/SpecDist1_new.pdf", width = 9, height = 5)
ggsave(compound_figure1, file = "results/SDM/SpecDist1_new.png", width = 9, height = 5)
ggsave(compound_figure2, file = "results/SDM/SpecDist2_new.pdf", width = 18, height = 10)
ggsave(compound_figure2, file = "results/SDM/SpecDist2_new.png", width = 18, height = 10)
'
