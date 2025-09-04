# This script downloads the weather data from Spartacus for the coordinates of the samples
# and concatenates the files into one file

# load the required libraries
library(dplyr)
library(lubridate)
library(terra)
library(tidyverse)
library(ggmap) # For creating maps with background tiles

# Set the working directory to the location of your script
setwd("D:/GitHub/UrbanDrosophilaEcology/")
# --- Set Working Directory ---
args <- commandArgs(trailingOnly = TRUE)
WD <- if (length(args) == 0) getwd() else args[1]
setwd(WD)
# setwd("/Users/martinkapun/Documents/GitHub/UrbanDrosophilaEcology")
# setwd("D:/GitHub/UrbanDrosophilaEcology")

dir.create("results/ClimateChange", showWarnings = FALSE)

# --- Data Import and Cleaning ---
DATA <- read.csv("data/Samples_inca_spartacus_vienna_clean_final.csv", header = TRUE)

# Only Viennese samples (omit rows with missing values)
DATA.Vienna <- na.omit(DATA)

# Combine the latitude and longitude into a single string
# and remove duplicates
#<host>/<version>/<type>/<mode>/<resource_id>
latlonglist <- DATA %>%
    select(Latitude, Longitude) %>%
    distinct() %>%
    mutate(lat_lon = paste0(Latitude, "%2C", Longitude)) %>%
    pull(lat_lon)
latlonglist <- unique(latlonglist)

# set the parameters for the API call
version <- "v1"
type <- "timeseries"
mode <- "historical"
resource_id <- "spartacus-v2-1y-1km"
parameters <- c("RR", "TM", "SA")
start <- "1990-01-01T00%3A00"
today <- Sys.Date()
end <- paste0(today, "T00%3A00")
host <- "https://dataset.api.hub.geosphere.at"

for (i in 1:length(latlonglist)) {
    LL <- paste0("lat_lon=", latlonglist[i])
    print(LL)

    # Create the URL
    url <- paste0(
        host, "/", version, "/", type, "/", mode, "/", resource_id,
        "?", paste(paste0("parameters=", parameters, collapse = "&")),
        "&start=", start,
        "&end=", end,
        "&", LL,
        # "&", paste(paste0("lat_lon=",latlonglist, collapse = "&")),
        "&output_format=csv"
    )

    # Download the data
    download.file(url, destfile = paste0("results/ClimateChange/Spartacus", i, ".csv"), mode = "wb", quiet = FALSE)
}

# Concatenate all the files into one file and remove the first row of each file,
# except the first one
file.create("results/ClimateChange/Spartacus_full.csv") # Create the file
for (i in seq_along(latlonglist)) {
    # Read the file
    temp <- read.csv(
        paste0("results/ClimateChange/Spartacus", i, ".csv"),
        header = TRUE,
        sep = ",",
        stringsAsFactors = FALSE
    )
    # Remove the first row of each file, except the first one
    if (i == 1) {
        write.csv(temp, file = "results/ClimateChange/Spartacus_full.csv", row.names = FALSE)
    } else {
        temp_file <- tempfile()
        write.csv(temp[-1, ], file = temp_file, row.names = FALSE)
        file.append("results/ClimateChange/Spartacus_full.csv", temp_file)
        unlink(temp_file)
    }
}
# Remove the individual files
for (i in seq_along(latlonglist)) {
    file.remove(paste0("results/ClimateChange/Spartacus", i, ".csv"))
}

## in the Spartacus_full.csv file, keep the header and remove all subsequent rows that have the same header
spartacus_full <- read.csv("results/ClimateChange/Spartacus_full.csv", header = TRUE, stringsAsFactors = FALSE)
spartacus_full <- spartacus_full[!duplicated(spartacus_full), ]

## for each unique combination of lat and lon, calculate generate a lineplot of the third column and use the year in the first colum as the x-axis
spartacus_full$Date <- as.Date(spartacus_full[, 1], format = "%Y-%m-%dT%H:%M+%S:00")
spartacus_full$Year <- year(spartacus_full$Date)
spartacus_full$Temp <- as.numeric(spartacus_full[, 3])

# Remove rows with missing temperature data
spartacus_full <- spartacus_full %>%
    filter(!is.na(Temp)) %>%
    filter(Year < 2025) # Exclude 2025 (incomplete year)
# 48.22669250879852, 16.328909082414857
# 48.172825544645214, 16.380336527494304
# Add Vienna district classification based on coordinates
# Inner Vienna districts (1-9) are roughly within these boundaries:
# Latitude: 48.195-48.225, Longitude: 16.35-16.395
spartacus_full <- spartacus_full %>%
    mutate(
        vienna_zone = case_when(
            lat >= 48.172825544645214 & lat <= 48.22669250879852 & lon >= 16.328909082414857 & lon <= 16.380336527494304 ~ "Inner_Vienna",
            TRUE ~ "Outer_Vienna"
        )
    )

# Check distribution
zone_summary <- spartacus_full %>%
    group_by(vienna_zone) %>%
    summarise(
        n_records = n(),
        n_locations = n_distinct(paste(lat, lon)),
        .groups = "drop"
    )

cat("Total valid temperature records:", nrow(spartacus_full), "\n")
cat("Vienna zone distribution:\n")
print(zone_summary)


## fit models to estimate temperature trends
# Linear trend model for overall data
model_linear <- lm(Temp ~ Year, data = spartacus_full)

# Zone-specific models
model_inner <- lm(Temp ~ Year, data = filter(spartacus_full, vienna_zone == "Inner_Vienna"))
model_outer <- lm(Temp ~ Year, data = filter(spartacus_full, vienna_zone == "Outer_Vienna"))

# Compare zone trends
inner_slope <- coef(model_inner)[2]
outer_slope <- coef(model_outer)[2]
inner_r2 <- summary(model_inner)$r.squared
outer_r2 <- summary(model_outer)$r.squared

cat("\n=== ZONE COMPARISON ===\n")
cat("Inner Vienna trend:", round(inner_slope, 4), "°C/year (R² =", round(inner_r2, 3), ")\n")
cat("Outer Vienna trend:", round(outer_slope, 4), "°C/year (R² =", round(outer_r2, 3), ")\n")
cat("Difference:", round(inner_slope - outer_slope, 4), "°C/year\n")

# Calculate temperature trends by location
# First check data availability per location
location_data_check <- spartacus_full %>%
    group_by(lat, lon, vienna_zone) %>%
    summarise(
        n_obs = n(),
        n_valid_temp = sum(!is.na(Temp)),
        n_years = n_distinct(Year),
        temp_range = ifelse(sum(!is.na(Temp)) > 1,
            max(Temp, na.rm = TRUE) - min(Temp, na.rm = TRUE),
            NA
        ),
        .groups = "drop"
    )

# Only analyze locations with sufficient data (at least 5 observations, 3 years, and valid temperature data)
locations_sufficient <- location_data_check %>%
    filter(n_obs >= 5, n_years >= 3, n_valid_temp >= 5, !is.na(temp_range))

cat(
    "Locations with sufficient data for trend analysis:", nrow(locations_sufficient),
    "out of", nrow(location_data_check), "\n"
)

# Calculate trends only for locations with sufficient data
if (nrow(locations_sufficient) > 0) {
    location_trends <- spartacus_full %>%
        filter(paste(lat, lon) %in% paste(locations_sufficient$lat, locations_sufficient$lon)) %>%
        group_by(lat, lon, vienna_zone) %>%
        summarise(
            n_obs = n(),
            slope = tryCatch(
                {
                    model <- lm(Temp ~ Year, data = cur_data())
                    coef(model)[2]
                },
                error = function(e) NA
            ),
            intercept = tryCatch(
                {
                    model <- lm(Temp ~ Year, data = cur_data())
                    coef(model)[1]
                },
                error = function(e) NA
            ),
            r_squared = tryCatch(
                {
                    model <- lm(Temp ~ Year, data = cur_data())
                    summary(model)$r.squared
                },
                error = function(e) NA
            ),
            p_value = tryCatch(
                {
                    model <- lm(Temp ~ Year, data = cur_data())
                    summary(model)$coefficients[2, 4]
                },
                error = function(e) NA
            ),
            .groups = "drop"
        ) %>%
        filter(!is.na(slope)) # Remove failed fits
} else {
    location_trends <- data.frame()
    cat("Warning: No locations have sufficient data for trend analysis
")
}

# Calculate temperature increase over the study period
year_range <- max(spartacus_full$Year) - min(spartacus_full$Year)

if (nrow(location_trends) > 0) {
    # Overall statistics
    temp_increase_stats <- location_trends %>%
        mutate(total_increase = slope * year_range) %>%
        summarise(
            avg_increase_per_year = mean(slope, na.rm = TRUE),
            min_increase_per_year = min(slope, na.rm = TRUE),
            max_increase_per_year = max(slope, na.rm = TRUE),
            avg_total_increase = mean(total_increase, na.rm = TRUE),
            min_total_increase = min(total_increase, na.rm = TRUE),
            max_total_increase = max(total_increase, na.rm = TRUE),
            avg_r2 = mean(r_squared, na.rm = TRUE),
            significant_trends = sum(p_value < 0.05, na.rm = TRUE),
            total_locations = n()
        )

    # Zone-specific statistics
    zone_stats <- location_trends %>%
        group_by(vienna_zone) %>%
        mutate(total_increase = slope * year_range) %>%
        summarise(
            n_locations = n(),
            avg_increase_per_year = mean(slope, na.rm = TRUE),
            avg_total_increase = mean(total_increase, na.rm = TRUE),
            avg_r2 = mean(r_squared, na.rm = TRUE),
            significant_trends = sum(p_value < 0.05, na.rm = TRUE),
            .groups = "drop"
        )
} else {
    # Use overall model if location-specific analysis fails
    overall_slope <- coef(model_linear)[2]
    temp_increase_stats <- data.frame(
        avg_increase_per_year = overall_slope,
        min_increase_per_year = overall_slope,
        max_increase_per_year = overall_slope,
        avg_total_increase = overall_slope * year_range,
        min_total_increase = overall_slope * year_range,
        max_total_increase = overall_slope * year_range,
        avg_r2 = summary(model_linear)$r.squared,
        significant_trends = ifelse(summary(model_linear)$coefficients[2, 4] < 0.05, 1, 0),
        total_locations = 1
    )

    # Create zone stats with overall model
    zone_stats <- data.frame(
        vienna_zone = c("Inner_Vienna", "Outer_Vienna"),
        n_locations = c(0, 1),
        avg_increase_per_year = c(inner_slope, outer_slope),
        avg_total_increase = c(inner_slope * year_range, outer_slope * year_range),
        avg_r2 = c(inner_r2, outer_r2),
        significant_trends = c(0, 0)
    )

    cat("Using overall linear model due to insufficient location-specific data\n")
}

# Print temperature trend statistics
cat("\n=== TEMPERATURE TREND ANALYSIS ===\n")
cat("Study period:", min(spartacus_full$Year), "-", max(spartacus_full$Year), "(", year_range, "years)\n")
cat("Average temperature increase per year:", round(temp_increase_stats$avg_increase_per_year, 4), "°C/year\n")
cat(
    "Range of temperature increase per year:", round(temp_increase_stats$min_increase_per_year, 4), "to",
    round(temp_increase_stats$max_increase_per_year, 4), "°C/year\n"
)
cat("Average total temperature increase:", round(temp_increase_stats$avg_total_increase, 2), "°C\n")
cat(
    "Range of total temperature increase:", round(temp_increase_stats$min_total_increase, 2), "to",
    round(temp_increase_stats$max_total_increase, 2), "°C\n"
)
cat("Average model R²:", round(temp_increase_stats$avg_r2, 3), "\n")
cat(
    "Locations with significant trends (p<0.05):", temp_increase_stats$significant_trends, "out of",
    temp_increase_stats$total_locations, "\n"
)

# Print zone-specific statistics
cat("\n=== ZONE-SPECIFIC TRENDS ===\n")
if (exists("zone_stats")) {
    for (i in 1:nrow(zone_stats)) {
        zone <- zone_stats$vienna_zone[i]
        cat(zone, ":\n")
        cat("  Locations:", zone_stats$n_locations[i], "\n")
        cat("  Avg increase:", round(zone_stats$avg_increase_per_year[i], 4), "°C/year\n")
        cat("  Total increase:", round(zone_stats$avg_total_increase[i], 2), "°C\n")
        cat("  Avg R²:", round(zone_stats$avg_r2[i], 3), "\n")
        cat("  Significant trends:", zone_stats$significant_trends[i], "\n\n")
    }
}

# Save trend statistics
if (nrow(location_trends) > 0) {
    write.csv(location_trends, "results/ClimateChange/temperature_trends_by_location.csv", row.names = FALSE)
}
write.csv(temp_increase_stats, "results/ClimateChange/temperature_increase_summary.csv", row.names = FALSE)
if (exists("zone_stats")) {
    write.csv(zone_stats, "results/ClimateChange/temperature_trends_by_zone.csv", row.names = FALSE)
}

PLOT <- spartacus_full %>%
    select(lat, lon, Year, Temp, vienna_zone) %>%
    rename(Latitude = lat, Longitude = lon) %>%
    arrange(Latitude, Longitude, Year) %>%
    ## line color by zone with transparency
    ggplot(aes(x = Year, y = Temp, group = interaction(Latitude, Longitude), color = vienna_zone)) +
    geom_line(alpha = 0.2, size = 0.3) +
    ## add zone-specific trend lines
    geom_smooth(aes(group = vienna_zone), method = "lm", se = TRUE, size = 1.2) +
    scale_color_manual(
        values = c("Inner_Vienna" = "red", "Outer_Vienna" = "blue"),
        name = "Vienna Zone",
        labels = if (exists("zone_stats")) {
            c(
                paste0(
                    "Inner Vienna (", round(zone_stats$avg_increase_per_year[zone_stats$vienna_zone == "Inner_Vienna"], 4),
                    " °C/yr, ", round(zone_stats$avg_total_increase[zone_stats$vienna_zone == "Inner_Vienna"], 2), " °C total)"
                ),
                paste0(
                    "Outer Vienna (", round(zone_stats$avg_increase_per_year[zone_stats$vienna_zone == "Outer_Vienna"], 4),
                    " °C/yr, ", round(zone_stats$avg_total_increase[zone_stats$vienna_zone == "Outer_Vienna"], 2), " °C total)"
                )
            )
        } else {
            c("Inner Vienna", "Outer Vienna")
        }
    ) +
    labs(
        x = "Year", y = "Temperature (°C)",
        title = "Temperature Trends by Vienna Zone",
        subtitle = paste0(
            "Overall: ", round(temp_increase_stats$avg_increase_per_year, 4),
            " °C/year (", round(temp_increase_stats$avg_total_increase, 2), " °C total)"
        )
    ) +
    theme_bw() +
    theme(legend.position = "bottom", legend.text = element_text(size = 10))

PLOT

# --- Vienna Map with Temperature Zones ---
# Create a summary of locations for mapping
locations_summary <- spartacus_full %>%
    group_by(lat, lon, vienna_zone) %>%
    summarise(
        n_records = n(),
        mean_temp = mean(Temp, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    rename(Latitude = lat, Longitude = lon)

# Define bounding box for Vienna map
bbox <- c(
    left = as.numeric(min(DATA.Vienna$Longitude)) - 0.02,
    bottom = as.numeric(min(DATA.Vienna$Latitude)) - 0.02,
    right = as.numeric(max(DATA.Vienna$Longitude)) + 0.02,
    top = as.numeric(max(DATA.Vienna$Latitude)) + 0.02
)

# Get Vienna map
vienna_map <- get_stadiamap(bbox = bbox, maptype = "stamen_toner", zoom = 11)

# Define inner Vienna coordinates for rectangle
inner_vienna_coords <- data.frame(
    lat_min = 48.172825544645214,
    lat_max = 48.22669250879852,
    lon_min = 16.328909082414857,
    lon_max = 16.380336527494304
)

# Create Vienna map with temperature monitoring locations
# Simple approach: use fixed size points to avoid scale issues
MAP_PLOT <- ggmap(vienna_map) +
    geom_point(
        data = locations_summary,
        aes(x = as.numeric(Longitude), y = as.numeric(Latitude), color = vienna_zone),
        size = 3, alpha = 0.8
    ) +
    # Add rectangle for inner Vienna zone
    geom_rect(
        aes(
            xmin = 16.328909082414857, xmax = 16.380336527494304,
            ymin = 48.172825544645214, ymax = 48.22669250879852
        ),
        fill = "red", alpha = 0.1, color = "red", linewidth = 1.5, inherit.aes = FALSE
    ) +
    scale_color_manual(
        values = c("Inner_Vienna" = "red", "Outer_Vienna" = "blue"),
        name = "Vienna Zone"
    ) +
    labs(
        title = "Vienna Temperature Monitoring Locations",
        subtitle = "Red rectangle shows Inner Vienna zone classification",
        x = "Longitude",
        y = "Latitude"
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom"
    )

MAP_PLOT

# Print model summary
cat("\n=== LINEAR TREND MODEL SUMMARY ===\n")
summary(model_linear)

# Save the plots
ggsave("results/ClimateChange/Temperature_Trends_by_Location.png", plot = PLOT, width = 10, height = 6, dpi = 300)
ggsave("results/ClimateChange/Temperature_Trends_by_Location.pdf", plot = PLOT, width = 10, height = 6, dpi = 300)
ggsave("results/ClimateChange/Vienna_Temperature_Map.png", plot = MAP_PLOT, width = 10, height = 8, dpi = 300)
ggsave("results/ClimateChange/Vienna_Temperature_Map.pdf", plot = MAP_PLOT, width = 10, height = 8, dpi = 300)

# Save the full dataset
