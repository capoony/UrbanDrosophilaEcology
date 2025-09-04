# Test script for Vienna map - minimal version to avoid data download
library(dplyr)
library(tidyverse)
library(ggmap)

# Set working directory
setwd("D:/GitHub/UrbanDrosophilaEcology/")

# Check if we have existing climate data
if (file.exists("results/ClimateChange/Spartacus_full.csv")) {
    cat("Loading existing climate data...\n")
    spartacus_full <- read.csv("results/ClimateChange/Spartacus_full.csv", header = TRUE, stringsAsFactors = FALSE)
    
    # Basic data processing
    spartacus_full$Date <- as.Date(spartacus_full[, 1], format = "%Y-%m-%dT%H:%M+%S:00")
    spartacus_full$Year <- lubridate::year(spartacus_full$Date)
    spartacus_full$Temp <- as.numeric(spartacus_full[, 3])
    
    # Remove rows with missing temperature data
    spartacus_full <- spartacus_full %>%
        filter(!is.na(Temp)) %>%
        filter(Year < 2025)
    
    # Add Vienna zone classification
    spartacus_full <- spartacus_full %>%
        mutate(
            vienna_zone = case_when(
                lat >= 48.172825544645214 & lat <= 48.22669250879852 & 
                lon >= 16.328909082414857 & lon <= 16.380336527494304 ~ "Inner_Vienna",
                TRUE ~ "Outer_Vienna"
            )
        )
    
    # Create locations summary
    locations_summary <- spartacus_full %>%
        group_by(lat, lon, vienna_zone) %>%
        summarise(
            n_records = n(),
            mean_temp = mean(Temp, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        rename(Latitude = lat, Longitude = lon)
    
    cat("Data summary:\n")
    print(head(locations_summary))
    cat("Range of n_records:", min(locations_summary$n_records), "to", max(locations_summary$n_records), "\n")
    cat("Unique vienna_zone values:", unique(locations_summary$vienna_zone), "\n")
    
} else {
    cat("No existing climate data found. Creating dummy data for testing...\n")
    # Create dummy data for testing
    locations_summary <- data.frame(
        Latitude = c(48.2, 48.21, 48.19, 48.15, 48.25),
        Longitude = c(16.37, 16.36, 16.38, 16.40, 16.35),
        vienna_zone = c("Inner_Vienna", "Inner_Vienna", "Inner_Vienna", "Outer_Vienna", "Outer_Vienna"),
        n_records = c(500, 300, 800, 200, 150),
        mean_temp = c(12.5, 12.3, 12.7, 12.1, 12.0)
    )
}

# Define bounding box for Vienna map
bbox <- c(
    left = min(locations_summary$Longitude) - 0.02,
    bottom = min(locations_summary$Latitude) - 0.02,
    right = max(locations_summary$Longitude) + 0.02,
    top = max(locations_summary$Latitude) + 0.02
)

cat("Bounding box:", bbox, "\n")

# Get Vienna map
tryCatch({
    vienna_map <- get_stadiamap(bbox = bbox, maptype = "stamen_toner", zoom = 10)
    cat("Map downloaded successfully.\n")
}, error = function(e) {
    cat("Error downloading map:", e$message, "\n")
    stop("Cannot proceed without map")
})

# Define inner Vienna coordinates for rectangle
inner_vienna_coords <- data.frame(
    lat_min = 48.172825544645214,
    lat_max = 48.22669250879852,
    lon_min = 16.328909082414857,
    lon_max = 16.380336527494304
)

# Create Vienna map with temperature monitoring locations
cat("Creating map plot...\n")

MAP_PLOT <- ggmap(vienna_map) +
    geom_point(
        data = locations_summary,
        aes(x = Longitude, y = Latitude, color = vienna_zone),
        size = 3, alpha = 0.8
    ) +
    # Add rectangle for inner Vienna zone
    geom_rect(
        data = inner_vienna_coords,
        aes(xmin = lon_min, xmax = lon_max, ymin = lat_min, ymax = lat_max),
        fill = "red", alpha = 0.2, color = "red", linewidth = 1.5, inherit.aes = FALSE
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

print(MAP_PLOT)

# Save the plot
dir.create("results/ClimateChange", showWarnings = FALSE)
ggsave("results/ClimateChange/Vienna_Temperature_Map_Test.png", plot = MAP_PLOT, width = 10, height = 8, dpi = 300)
cat("Map saved successfully!\n")
