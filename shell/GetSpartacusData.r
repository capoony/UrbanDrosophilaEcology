library(ncdf4)
library(terra)
library(httr)
library(dplyr)
library(lubridate)

## bounding box
vmaxLat <- 48.5
vminLat <- 47.5
vminLon <- 15.5
vmaxLon <- 16.7

# --- Set Working Directory ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  WD <- getwd()
} else {
  WD <- args[1]
}
setwd(WD)

# Base URL for Geosphere Austria
base_url <- "https://public.hub.geosphere.at/datahub/resources/spartacus-v2-"

# Bounding box (WGS84 coordinates)
bbox <- list(min_lon = vminLon, min_lat = vminLat, max_lon = vmaxLon, max_lat = vmaxLat)

# Function to download a file
download_file <- function(url, save_path) {
  tryCatch(
    {
      GET(url, write_disk(save_path, overwrite = TRUE))
      message("Downloaded: ", save_path)
    },
    error = function(e) {
      message("Failed to download ", url, ": ", e$message)
    }
  )
}

# Function to reproject, crop, interpolate, and print extent of a NetCDF file
process_netcdf <- function(file_path) {
  tryCatch(
    {
      # Open the NetCDF file as a SpatRaster
      r <- terra::rast(file_path)

      # Reproject to WGS84
      r <- terra::project(r, "+proj=longlat +datum=WGS84")

      # Print extent before cropping
      ext <- ext(r)
      message("Reprojected file extent: min_lon=", ext[1], ", min_lat=", ext[3], ", max_lon=", ext[2], ", max_lat=", ext[4])

      # Crop to bounding box
      r <- crop(r, ext(bbox$min_lon, bbox$max_lon, bbox$min_lat, bbox$max_lat))

      # Interpolate to 100x100m grid cells (approx. 0.0009 degrees)
      r <- resample(r, rast(ext(r), res = c(0.0009, 0.0009), crs = "+proj=longlat +datum=WGS84"), method = "bilinear")

      # Print extent after cropping and interpolation
      ext <- ext(r)
      message("Processed file extent: min_lon=", ext[1], ", min_lat=", ext[3], ", max_lon=", ext[2], ", max_lat=", ext[4])

      # Save the processed file
      processed_path <- gsub("\\.nc$", "_processed.tif", file_path)
      writeRaster(r, processed_path, overwrite = TRUE)
      message("Processed and saved: ", processed_path)
    },
    error = function(e) {
      message("Failed to process ", file_path, ": ", e$message)
    }
  )
}

# Parameters
output_folder <- "data/spartacus_data"
year <- 2024
params <- list(
  yearly = list(suffix = "YEARLY_", variables = c("TM", "RR", "SA")),
  monthly = list(suffix = "MONTHLY_", variables = c("TM", "RR", "SA")),
  daily = list(suffix = "DAILY_", variables = c("TX", "TN", "RR", "SA"))
)

# Create output directory
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# Loop through each time period and parameter
for (period in names(params)) {
  period_dir <- file.path(output_folder, period)
  dir.create(period_dir, showWarnings = FALSE)

  for (var in params[[period]]$variables) {
    filename <- sprintf("SPARTACUS2-%s%s_%d.nc", params[[period]]$suffix, var, year)
    period_suffix <- switch(period,
      yearly = "1y-",
      monthly = "1m-",
      daily = "1d-"
    )
    url <- sprintf("%s%s1km/filelisting/%s/%s", base_url, period_suffix, var, filename)
    save_path <- file.path(period_dir, filename)

    # Download file
    download_file(url, save_path)

    # Process file
    process_netcdf(save_path)
  }
}

message("All downloads and processing completed.")

#### now extract values for coordinates and Average 14days prior to collection date for daily data

# Parameters
processed_tif_directory <- "data/spartacus_data" # Directory with processed TIFFs
input_csv <- "data/Samples_inca.csv" # CSV file with coordinates and collectionEnd
tif_patterns <- list(
  yearly = "SPARTACUS2-YEARLY_.*_processed.tif$",
  monthly = "SPARTACUS2-MONTHLY_.*_processed.tif$",
  daily = "SPARTACUS2-DAILY_.*_processed.tif$"
)
output_csv <- "data/Samples_inca_spartacus.csv"

# Read CSV
coordinates <- read.csv(input_csv, header = TRUE)
coordinates$collectionEnd <- dmy(coordinates$collectionEnd)

# Function to extract value for yearly data
extract_yearly <- function(r, coord) {
  unlist(extract(r, cbind(coord$Longitude, coord$Latitude)))[1]
}

# Function to extract value for monthly data (returns only the sampling month value)
extract_monthly <- function(r, coord, collection_date) {
  sampling_month <- month(collection_date)
  r_time <- time(r)
  idx <- which(month(r_time) == sampling_month)
  if (length(idx) == 0) {
    return(NA)
  }
  values <- extract(r[[idx]], cbind(coord$Longitude, coord$Latitude))
  if (all(is.na(values))) {
    return(NA)
  }
  return(unlist(values)[1])
}

# Function to extract value for daily data
extract_daily <- function(r, coord, collection_date) {
  if (is.na(collection_date)) {
    return(NA)
  } # Handle missing dates
  sampling_range <- seq(as.Date(collection_date) - 14, as.Date(collection_date), by = "1 day")
  r_time <- time(r)
  if (is.null(r_time)) {
    warning("Raster does not have valid time metadata.")
    return(NA)
  }
  idx <- which(as.Date(r_time) %in% sampling_range)
  if (length(idx) == 0) {
    warning("No matching dates found in the raster for the sampling range.")
    return(NA)
  }
  values <- extract(r[[idx]], cbind(coord$Longitude, coord$Latitude))
  if (all(is.na(values))) {
    warning("All extracted values are NA for the given range.")
    return(NA)
  }
  return(mean(unlist(values), na.rm = TRUE))
}

# Extract values
DATA <- coordinates %>%
  rowwise() %>%
  mutate(
    Yearly_TM = {
      tif_files <- list.files(file.path(processed_tif_directory, "yearly"), pattern = "SPARTACUS2-YEARLY_TM_.*_processed.tif$", full.names = TRUE)
      if (length(tif_files) == 0) {
        return(NA)
      }
      r <- rast(tif_files[1])
      extract_yearly(r, cur_data())
    },
    Yearly_RR = {
      tif_files <- list.files(file.path(processed_tif_directory, "yearly"), pattern = "SPARTACUS2-YEARLY_RR_.*_processed.tif$", full.names = TRUE)
      if (length(tif_files) == 0) {
        return(NA)
      }
      r <- rast(tif_files[1])
      extract_yearly(r, cur_data())
    },
    Yearly_SA = {
      tif_files <- list.files(file.path(processed_tif_directory, "yearly"), pattern = "SPARTACUS2-YEARLY_SA_.*_processed.tif$", full.names = TRUE)
      if (length(tif_files) == 0) {
        return(NA)
      }
      r <- rast(tif_files[1])
      extract_yearly(r, cur_data())
    },
    Monthly_TM = {
      tif_files <- list.files(file.path(processed_tif_directory, "monthly"), pattern = "SPARTACUS2-MONTHLY_TM_.*_processed.tif$", full.names = TRUE)
      if (length(tif_files) == 0) {
        return(NA)
      }
      r <- rast(tif_files[1])
      extract_monthly(r, cur_data(), collectionEnd)
    },
    Monthly_RR = {
      tif_files <- list.files(file.path(processed_tif_directory, "monthly"), pattern = "SPARTACUS2-MONTHLY_RR_.*_processed.tif$", full.names = TRUE)
      if (length(tif_files) == 0) {
        return(NA)
      }
      r <- rast(tif_files[1])
      extract_monthly(r, cur_data(), collectionEnd)
    },
    Monthly_SA = {
      tif_files <- list.files(file.path(processed_tif_directory, "monthly"), pattern = "SPARTACUS2-MONTHLY_SA_.*_processed.tif$", full.names = TRUE)
      if (length(tif_files) == 0) {
        return(NA)
      }
      r <- rast(tif_files[1])
      extract_monthly(r, cur_data(), collectionEnd)
    },
    Daily_TX = {
      tif_files <- list.files(file.path(processed_tif_directory, "daily"), pattern = "SPARTACUS2-DAILY_TX_.*_processed.tif$", full.names = TRUE)
      if (length(tif_files) == 0) {
        return(NA)
      }
      r <- rast(tif_files[1])
      extract_daily(r, cur_data(), collectionEnd)
    },
    Daily_TN = {
      tif_files <- list.files(file.path(processed_tif_directory, "daily"), pattern = "SPARTACUS2-DAILY_TN_.*_processed.tif$", full.names = TRUE)
      if (length(tif_files) == 0) {
        return(NA)
      }
      r <- rast(tif_files[1])
      extract_daily(r, cur_data(), collectionEnd)
    },
    Daily_RR = {
      tif_files <- list.files(file.path(processed_tif_directory, "daily"), pattern = "SPARTACUS2-DAILY_RR_.*_processed.tif$", full.names = TRUE)
      if (length(tif_files) == 0) {
        return(NA)
      }
      r <- rast(tif_files[1])
      extract_daily(r, cur_data(), collectionEnd)
    },
    Daily_SA = {
      tif_files <- list.files(file.path(processed_tif_directory, "daily"), pattern = "SPARTACUS2-DAILY_SA_.*_processed.tif$", full.names = TRUE)
      if (length(tif_files) == 0) {
        return(NA)
      }
      r <- rast(tif_files[1])
      extract_daily(r, cur_data(), collectionEnd)
    }
  ) %>%
  ungroup()

# Write to CSV
write.table(as.matrix(DATA), output_csv, row.names = FALSE, quote = FALSE, sep = ",")
message("Extracted values saved to: ", output_csv)
