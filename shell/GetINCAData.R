# Install packages if needed:
# install.packages("terra")

library(terra)

# --- Set Working Directory ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    WD <- getwd()
} else {
    WD <- args[1]
}
setwd(WD)

# List of variables (the part that replaces '*'):
varList <- c("GL", "RH2M", "RR", "T2M", "TD2M", "UU", "VV", "P0")

# GL: Globalstrahlung (Global Solar Radiation) – gemessene oder berechnete Gesamtmenge an Sonnenstrahlung, die auf eine horizontale Fläche trifft (Einheit: W/m² oder kWh/m²).

# RH2M: Relative Luftfeuchtigkeit in 2 Metern Höhe (Relative Humidity at 2 Meters) – prozentualer Anteil der tatsächlichen Wasserdampfmenge zur maximal möglichen in der Luft (Einheit: %).

# RR: Niederschlagssumme (Rainfall Rate) – gesamte Niederschlagsmenge in einem bestimmten Zeitraum (Einheit: mm oder kg/m²).

# T2M: Lufttemperatur in 2 Metern Höhe (Temperature at 2 Meters) – gemessene oder berechnete Temperatur in einer Standardhöhe von 2 Metern über dem Boden (Einheit: °C).

# TD2M: Taupunkttemperatur in 2 Metern Höhe (Dew Point Temperature at 2 Meters) – Temperatur, bei der die Luftfeuchtigkeit 100 % erreicht und Kondensation einsetzt (Einheit: °C).

# UU: Windgeschwindigkeit (Wind Speed) – gemessene oder berechnete Geschwindigkeit des Windes (Einheit: m/s).

# VV: Windrichtung (Wind Direction) – gemessene oder berechnete Richtung, aus der der Wind kommt, angegeben in Grad (0° = Norden, 90° = Osten, 180° = Süden, 270° = Westen).

# Months from June to December 2024 (202406, 202407, ..., 202412)
monthList <- sprintf("2024%02d", 6:12)

# Create folders to store downloaded files and outputs
dir.create("data/downloads", showWarnings = FALSE)
dir.create("data/outputs_full", showWarnings = FALSE)
dir.create("data/outputs", showWarnings = FALSE)

###############################################################################
# Helper function to concatenate a list of SpatRasters in a loop
###############################################################################
concatSpatRasters <- function(raster_list) {
    # If there's only one element, just return it
    if (length(raster_list) == 1) {
        return(raster_list[[1]])
    }
    # Otherwise, start with the first
    result <- raster_list[[1]]
    # Concatenate subsequent rasters
    for (i in seq(2, length(raster_list))) {
        result <- c(result, raster_list[[i]])
    }
    return(result)
}

# Loop over each variable
for (v in varList) {
    # We'll collect the daily-averaged layers for all months in this list
    all_months_list <- list()

    # Loop over each month
    for (m in monthList) {
        # Construct the download URL, e.g.:
        # https://public.hub.geosphere.at/datahub/resources/inca-v1-1h-1km/filelisting/VV/INCAL_HOURLY_VV_202406.nc
        url <- paste0(
            "https://public.hub.geosphere.at/datahub/resources/inca-v1-1h-1km/filelisting/",
            v, "/INCAL_HOURLY_", v, "_", m, ".nc"
        )

        # Local filename for the downloaded file
        localfile <- file.path("data/downloads", paste0("INCAL_HOURLY_", v, "_", m, ".nc"))

        # Download the file (mode = "wb" ensures proper binary writing on Windows)
        cat("Downloading", url, "...\n")
        download.file(url, destfile = localfile, mode = "wb", quiet = FALSE)

        # Open the downloaded netCDF as a SpatRaster
        cat("Reading:", localfile, "\n")
        r <- rast(localfile) # Each layer corresponds to one hour

        # Define a daily chunk size in hours
        chunk_size <- 24 # 24 hours = 1 day
        n_hours <- nlyr(r) # total number of hourly layers in the file

        # Identify the start indices for each day
        chunk_starts <- seq(1, n_hours, by = chunk_size)

        # We'll store the daily mean for each day in a list
        monthly_means_list <- list()

        # Create a progress bar for daily chunks
        pb <- txtProgressBar(min = 0, max = length(chunk_starts), style = 3)

        for (i in seq_along(chunk_starts)) {
            start_idx <- chunk_starts[i]
            end_idx <- min(start_idx + chunk_size - 1, n_hours)

            # Subset the hourly layers for this day
            sub_r <- r[[start_idx:end_idx]]

            # Average across those hours (produces a single-layer SpatRaster)
            day_mean <- mean(sub_r)

            # Name the layer to reflect the day
            layer_name <- paste0(m, "_hours_", start_idx, "_to_", end_idx)
            names(day_mean) <- layer_name

            # Store in the monthly means list
            monthly_means_list[[layer_name]] <- day_mean

            # Update the progress bar
            setTxtProgressBar(pb, i)
        }

        # Close the progress bar
        close(pb)

        # Combine all daily means for this month
        # Concatenate this month's daily rasters into a single multi-layer SpatRaster
        if (length(monthly_means_list) > 0) {
            monthly_means_stack <- concatSpatRasters(monthly_means_list)
        } else {
            # If no layers, skip or create an empty SpatRaster
            monthly_means_stack <- NULL
        }

        # Collect for final combination (across all months)
        all_months_list[[m]] <- monthly_means_stack

        # Delete the downloaded file to save space
        cat("Removing file:", localfile, "\n")
        file.remove(localfile)
    }

    # Combine all months’ daily means for the current variable
    all_months_list <- Filter(Negate(is.null), all_months_list)

    if (length(all_months_list) > 0) {
        final_stack <- concatSpatRasters(all_months_list)
    } else {
        # No data at all
        final_stack <- NULL
    }



    # If final_stack is not NULL, proceed to write
    if (!is.null(final_stack)) {
        # Construct output file name
        out_file <- file.path("data/outputs", paste0("INCAL_dailyMean_", v, "_202406_202412.nc"))

        cat("Writing final daily means for variable", v, "to:", out_file, "\n")
        writeCDF(final_stack, filename = out_file, overwrite = TRUE)
    } else {
        cat("No data found for variable", v, "from June–December 2024.\n")
    }
}

cat("Done!\n")

### now restrict datasets to Vienna, interpolate to 100x100m
library(ncdf4)
library(raster)
library(rgdal)
library(sp)
library(gstat)

custom_tempdir <- "/path/to/custom/tempdir"

# Create the directory if it doesn't exist
if (!dir.exists(custom_tempdir)) {
    dir.create(custom_tempdir, recursive = TRUE)
}

# Set the custom directory as the R temporary directory
Sys.setenv(TMPDIR = custom_tempdir)


# Define bounding box
vmaxLat <- 48.5
vminLat <- 47.5
vminLon <- 15.5
vmaxLon <- 16.7
bbox <- extent(vminLon, vmaxLon, vminLat, vmaxLat)

# Folder containing NetCDF files
input_folder <- "/media/inter/mkapun/projects/RDA_VCF/data/outputs/"
output_folder <- "/media/inter/mkapun/projects/RDA_VCF/data/outputs_full/"

# Get list of all NetCDF files starting with "INCAL"
nc_files <- list.files(input_folder, pattern = "^INCAL.*\\.nc$", full.names = TRUE)

# Loop through each file
for (nc_file in nc_files) {
    # Load the NetCDF file
    nc_data <- brick(nc_file)

    # Reproject to WGS 84
    reprojected_data <- projectRaster(nc_data, crs = CRS("+proj=longlat +datum=WGS84"))

    # Crop to the bounding box
    cropped_data <- crop(reprojected_data, bbox)

    # Create a 100x100m grid
    res <- 0.001 # ~100m resolution in degrees
    target_grid <- raster(extent(bbox), crs = CRS("+proj=longlat +datum=WGS84"))
    res(target_grid) <- res

    # Interpolate to the finer grid
    interpolated_data <- resample(cropped_data, target_grid, method = "bilinear")

    # Rename layers to consecutive daily dates starting from 20240601
    start_date <- as.Date("2024-06-01")
    num_layers <- nlayers(interpolated_data)
    layer_names <- as.character(seq(start_date, by = "day", length.out = num_layers))
    names(interpolated_data) <- layer_names

    # Create the output filename
    base_name <- gsub("\\.nc$", "", basename(nc_file))
    output_file <- file.path(output_folder, paste0(base_name, "_1d_.nc"))

    # Save the output
    writeRaster(interpolated_data, output_file, format = "CDF", overwrite = TRUE)

    cat("Processed and saved:", output_file, "\n")
}
