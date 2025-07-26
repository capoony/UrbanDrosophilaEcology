# ------------------------------------------------------------------------------
# Script: GetViennaDataCubeData.r
# Purpose:
#   - Copy climate GeoTIFF files to a local directory
#   - Reproject all TIFF files to EPSG:4326
#   - Extract raster values at sample coordinates from a CSV
#   - Merge extracted values with the CSV and save the result
# Author: Martin Kapun
# ------------------------------------------------------------------------------

# --- Step 0: Copy Climate Data Files ---
# (Assumes this is run in a shell, not R)
# cp ~/mounts/BioMem_2/ssteindl/UC3/ClimateData/Vienna/S3/ViennaData/fairicube/100m/*/*.tif data/Climate

# --- Step 1: Load Required Libraries ---
library(raster) # For raster data handling
library(terra) # For raster operations (modern alternative to raster)
library(sf) # For spatial vector data
library(dplyr) # For data manipulation
library(tidyr) # For data reshaping
library(tools) # For file path manipulation

# --- Step 2: Set Working Directory ---
args <- commandArgs(trailingOnly = TRUE)
WD <- if (length(args) == 0) getwd() else args[1]
setwd(WD)

# --- Step 3: List and Prepare TIFF Files ---
tiff_dir <- "data/Climate/"
tiff_files <- list.files(path = tiff_dir, pattern = "\\.tif$", full.names = TRUE)

# Define the target CRS (EPSG:4326)
crs_epsg_4326 <- "EPSG:4326"

# Create output directory for reprojected files
output_dir <- "data/Climate/Reprojected/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Step 4: Reproject TIFF Files to EPSG:4326 ---
lapply(tiff_files, function(file) {
    tryCatch(
        {
            cat("\nProcessing file:", file, "\n")
            raster_file <- terra::rast(file)
            # If CRS is missing, assign the source CRS manually
            if (is.na(terra::crs(raster_file))) {
                cat("CRS is missing for file:", file, "- Assigning custom CRS\n")
                terra::crs(raster_file) <- "+proj=tmerc +lat_0=0 +lon_0=16.3333333333333 +k=1 +x_0=0 +y_0=-5000000 +ellps=bessel +towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.42319999999019 +units=m +no_defs"
            }
            # Reproject raster
            raster_reprojected <- terra::project(raster_file, crs_epsg_4326)
            # Save reprojected raster
            output_file <- file.path(output_dir, paste0("reprojected_", basename(file)))
            terra::writeRaster(raster_reprojected, output_file, overwrite = TRUE)
            cat("Saved reprojected file to:", output_file, "\n")
            return(output_file)
        },
        error = function(e) {
            cat("Error processing file:", file, "\nError message:", e$message, "\n")
            return(NULL)
        }
    )
})

# --- Step 5: Extract Raster Values at Sample Points ---

# Define file paths
geotiff_folder <- "data/Climate/Reprojected"
csv_path <- "data/Samples_inca_spartacus.csv"
output_csv_path <- "data/Samples_inca_spartacus_vienna.csv"

# Load sample data CSV
csv_data <- read.csv(csv_path)

# Extract coordinates (assumes Latitude is column 5, Longitude is column 6)
coords <- csv_data %>%
    select(Latitude = 5, Longitude = 6)

# Convert coordinates to sf object (WGS84)
points_sf <- st_as_sf(coords, coords = c("Longitude", "Latitude"), crs = 4326)

# List all reprojected GeoTIFF files
tif_files <- list.files(geotiff_folder, pattern = "\\.tif$", full.names = TRUE)

# Initialize list to store extracted values
extracted_values <- list()

# Loop through each GeoTIFF and extract values at sample points
for (tif in tif_files) {
    raster_stack <- raster::stack(tif)
    values <- raster::extract(raster_stack, st_coordinates(points_sf), method = "simple")
    values_df <- as.data.frame(values)
    colnames(values_df) <- paste0("Band_", seq_along(values_df))
    file_name <- tools::file_path_sans_ext(basename(tif))
    colnames(values_df) <- paste(file_name, colnames(values_df), sep = "_")
    extracted_values[[tif]] <- values_df
}

# Combine all extracted values into a single data frame
all_extracted_values <- bind_cols(extracted_values)

# --- Step 6: Rename Columns for Land Use Bands (if applicable) ---

# Define land use categories (for a specific raster file)
categories <- c(
    "looseresidential", "mediumresidential", "denseresidential", "largeisolatedresidential",
    "office", "commercial", "businesscore", "lowmixeduse", "industrial", "leisure", "health",
    "education", "indoorsports", "military", "sewage", "energy", "water", "construction",
    "greenstreets", "non-greenstreets", "parking", "railway", "transport", "park",
    "outdoorsports", "cemetery", "fields", "vineyards", "nursery", "forest", "meadow", "water"
)

# Update column names for the specific raster file
all_extracted_values <- all_extracted_values %>%
    rename_with(
        .fn = function(col) {
            pattern <- "reprojected_r01_real_land_use2020_100m_b32_1_1_Band_(\\d+)$"
            if (grepl(pattern, col)) {
                band_number <- as.numeric(sub(pattern, "\\1", col))
                if (!is.na(band_number) && band_number <= length(categories)) {
                    sub(
                        paste0("_Band_", band_number, "$"),
                        paste0("_", categories[band_number]),
                        col
                    )
                } else {
                    col
                }
            } else {
                col
            }
        }
    )

# Print updated column names for verification
print(colnames(all_extracted_values))

# --- Step 7: Merge Extracted Values with Original CSV and Save ---

updated_csv_data <- bind_cols(csv_data, all_extracted_values)
write.csv(updated_csv_data, output_csv_path, row.names = FALSE, quote = FALSE)
cat("Updated CSV saved to:", output_csv_path, "\n")
