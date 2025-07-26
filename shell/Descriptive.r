# -------------------------------------------------------------------
# Descriptive Analysis of Urban Drosophila Ecology Data
#
# This script performs basic data exploration, visualization, and
# summary statistics for Drosophila samples collected in Vienna and
# other locations. It generates:
#   - Sample counts per participant and per month
#   - Maps of sampling locations
#   - Species abundance barplots
#   - Temporal histograms of species counts
#
# Usage:
#   Rscript Descriptive.r [working_directory]
#
# Input:
#   data/Samples_inca_spartacus_vienna_clean.csv
#
# Output:
#   results/Descriptive/*.pdf (various plots)
#
# Author: Martin Kapun
# -------------------------------------------------------------------

# --- Load Required Libraries ---
library(tidyverse) # Data manipulation (dplyr, mutate, group_by, etc.) and visualization (ggplot2)
library(ggmap) # Fetching and plotting background maps (get_stadiamap, ggmap)
library(patchwork) # Combining multiple ggplot2 plots (plot_layout, / operator)
library(osmdata) # For getbb() (bounding box from OpenStreetMap)
library(FactoMineR) # Multivariate analysis (not used here, but often for PCA)
library(factoextra) # Visualization for FactoMineR outputs
library(vegan) # Ecological analysis (RDA, diversity, etc.)
library(terra) # Raster and spatial data (TIFFs)
library(raster) # Raster data handling
library(stars) # Spatiotemporal arrays
library(sf) # Simple features for spatial data
library(ggnewscale) # Multiple fill/color scales in ggplot2
library(ggpubr) # Publication-ready ggplot2 visualizations
library(lubridate) # Date/time manipulation
library(RColorBrewer) # Color palettes

# --- Set Working Directory ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  WD <- getwd()
} else {
  WD <- args[1]
}
setwd(WD)

# --- Data Loading and Preparation ---
# Read cleaned sample data
DATA <- read.csv("data/Samples_inca_spartacus_vienna_clean.csv", header = TRUE)

# Remove non-Viennese factors and omit rows with missing values
DATA.all <- DATA[, 1:36] %>% na.omit()

# Only Viennese samples (omit rows with missing values)
DATA.Vienna <- na.omit(DATA)

# Create output directory for results
dir.create("results/Descriptive", showWarnings = FALSE)

# --- Participant Summary ---
# Summarize number of samples per participant
sum_part <- DATA.Vienna %>%
  group_by(ParticipantId) %>%
  summarize(Count = n()) %>%
  arrange(desc(Count)) %>%
  group_by(Count) %>%
  summarize(Number = n()) %>%
  arrange(desc(Number))

# --- Month Translation (German to English) ---
month_translation <- c(
  "Januar" = "January", "Februar" = "February", "März" = "March",
  "April" = "April", "Mai" = "May", "Juni" = "June", "Juli" = "July",
  "August" = "August", "September" = "September", "Oktober" = "October",
  "November" = "November", "Dezember" = "December"
)

# --- Plot: Samples per Month (Vienna) ---
DATA.sub <- DATA.Vienna %>%
  mutate(
    collectionEnd = as.Date(collectionEnd, format = "%d/%m/%Y"),
    german_month_name = format(collectionEnd, "%B"),
    month_name = month_translation[german_month_name]
  )
samples_per_month <- DATA.sub %>%
  count(month_name) %>%
  mutate(month_name = factor(month_name, levels = month.name))
PLOT <- ggplot(samples_per_month, aes(x = month_name, y = n)) +
  geom_bar(stat = "identity", fill = "grey", color = "black") +
  theme_bw() +
  labs(
    title = "Number of Samples per Month",
    x = "Month",
    y = "Number of Samples"
  )
ggsave("results/Descriptive/SamplesPerMonth.Vienna.pdf", PLOT, width = 10, height = 8)

# --- Plot: Samples per Month (All) ---
DATA.sub <- DATA.all %>%
  mutate(
    collectionEnd = as.Date(collectionEnd, format = "%d/%m/%Y"),
    german_month_name = format(collectionEnd, "%B"),
    month_name = month_translation[german_month_name]
  )
samples_per_month <- DATA.sub %>%
  count(month_name) %>%
  mutate(month_name = factor(month_name, levels = month.name))
PLOT <- ggplot(samples_per_month, aes(x = month_name, y = n)) +
  geom_bar(stat = "identity", fill = "grey", color = "black") +
  theme_bw() +
  labs(
    title = "Number of Samples per Month",
    x = "Month",
    y = "Number of Samples"
  )
ggsave("results/Descriptive/SamplesPerMonth.all.pdf", PLOT, width = 10, height = 8)

# --- Map: Sampling Locations (All) ---
data_summary <- DATA.all[, 5:6] %>%
  group_by(Latitude, Longitude) %>%
  summarise(Samples = n(), .groups = "drop") %>%
  mutate(Sample_Class = cut(
    Samples,
    breaks = c(0, 1, 5, Inf), labels = c("1", "1-5", ">5")
  ))
size_mapping <- c("1" = 1, "1-5" = 2, ">5" = 3)
data_summary$Dot_Size <- size_mapping[as.character(data_summary$Sample_Class)]
bbox <- c(
  left = min(data_summary$Longitude) - 0.05,
  bottom = min(data_summary$Latitude) - 0.05,
  right = max(data_summary$Longitude) + 0.05,
  top = max(data_summary$Latitude) + 0.05
)
map <- get_stadiamap(bbox = bbox, maptype = "stamen_toner", zoom = 11)
PLOT <- ggmap(map) +
  geom_point(
    data = data_summary, aes(x = Longitude, y = Latitude, size = Sample_Class),
    color = "blue", alpha = 0.7
  ) +
  scale_size_manual(
    name = "Sample Class",
    values = c("1" = 1, "1-5" = 2, ">5" = 3),
    labels = c("1 Sample", "1-5 Samples", ">5 Samples")
  ) +
  theme_minimal() +
  labs(
    title = "Sampling Locations with Sample Class",
    x = "Longitude", y = "Latitude"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()
ggsave("results/Descriptive/Map_all.pdf", PLOT, width = 10, height = 8)

# --- Map: Sampling Locations (Vienna) ---
data_summary <- DATA.Vienna[, 5:6] %>%
  group_by(Latitude, Longitude) %>%
  summarise(Samples = n(), .groups = "drop") %>%
  mutate(Sample_Class = cut(
    Samples,
    breaks = c(0, 1, 5, Inf), labels = c("1", "1-5", ">5")
  ))
size_mapping <- c("1" = 2, "1-5" = 3, ">5" = 4)
data_summary$Dot_Size <- size_mapping[as.character(data_summary$Sample_Class)]
bbox <- c(
  left = min(data_summary$Longitude) - 0.05,
  bottom = min(data_summary$Latitude) - 0.05,
  right = max(data_summary$Longitude) + 0.05,
  top = max(data_summary$Latitude) + 0.05
)
map <- get_stadiamap(bbox = bbox, maptype = "stamen_toner", zoom = 11)
PLOT <- ggmap(map) +
  geom_point(
    data = data_summary, aes(x = Longitude, y = Latitude, size = Sample_Class),
    color = "blue", alpha = 0.7
  ) +
  scale_size_manual(
    name = "Sample Class",
    values = c("1" = 1, "1-5" = 3, ">5" = 5),
    labels = c("1 Sample", "1-5 Samples", ">5 Samples")
  ) +
  theme_minimal() +
  labs(
    title = "Sampling Locations with Sample Class",
    x = "Longitude", y = "Latitude"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()
ggsave("results/Descriptive/Map_Vienna.pdf", PLOT, width = 10, height = 8)

# --- Define Color Palette for Species ---
wes_distinct_palette <- c(
  "#EAD2AC", "#C29E73", "#8F5A5A", "#5E8D88", "#2F4B4F",
  "#F3B061", "#F05E23", "#D9A5A5", "#87AFC7", "#486F8A",
  "#A8BDBE", "#7C878F", "#D4AF37"
)

# --- Species Abundance Barplots (Vienna) ---
species_data <- DATA.Vienna[, 7:19]
species_names <- colnames(species_data)
species_colors <- setNames(wes_distinct_palette[1:length(species_names)], species_names)
species_totals <- colSums(species_data, na.rm = TRUE)
species_means <- apply(species_data, 2, function(x) mean(x[x > 0], na.rm = TRUE))
species_sd <- apply(species_data, 2, function(x) sd(x[x > 0], na.rm = TRUE))
species_n <- apply(species_data, 2, function(x) sum(x > 0))
species_se <- species_sd / sqrt(species_n)
species_df <- data.frame(
  Species = colnames(species_data),
  TotalAbundance = species_totals,
  MeanAbundance = species_means,
  SE = species_se
) %>%
  arrange(desc(TotalAbundance)) %>%
  mutate(Species = factor(Species, levels = Species)) %>%
  mutate(MeanAbundance = ifelse(MeanAbundance <= 0, 0.01, MeanAbundance)) %>%
  mutate(
    FlippedMeanAbundance = MeanAbundance,
    FlippedSEMin = MeanAbundance - SE,
    FlippedSEMax = MeanAbundance + SE
  )
total_abundance_plot <- ggplot(species_df, aes(x = Species, y = TotalAbundance, fill = Species)) +
  geom_bar(stat = "identity") +
  scale_y_log10() +
  scale_fill_manual(values = species_colors) +
  theme_bw() +
  labs(y = "Log Total Abundance", x = NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  guides(fill = "none")
relative_abundance_plot <- ggplot(species_df, aes(x = Species, y = MeanAbundance, fill = Species)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(
    ymin = pmax(FlippedSEMin, 0.01),
    ymax = FlippedSEMax
  ), width = 0.2) +
  scale_fill_manual(values = species_colors) +
  scale_y_reverse(breaks = seq(0, max(species_df$MeanAbundance, na.rm = TRUE), by = 5)) +
  theme_bw() +
  labs(y = "Mean Abundance (Non-Zero, with SE)", x = NULL) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  guides(fill = "none")
combined_plot <- total_abundance_plot / relative_abundance_plot +
  plot_layout(heights = c(1, 3))
ggsave("results/Descriptive/Abundance.Vienna.pdf", combined_plot, width = 10, height = 10)

# --- Species Abundance Barplots (All) ---
species_data <- DATA.all[, 7:19]
species_names <- colnames(species_data)
species_colors <- setNames(wes_distinct_palette[1:length(species_names)], species_names)
species_totals <- colSums(species_data, na.rm = TRUE)
species_means <- apply(species_data, 2, function(x) mean(x[x > 0], na.rm = TRUE))
species_sd <- apply(species_data, 2, function(x) sd(x[x > 0], na.rm = TRUE))
species_n <- apply(species_data, 2, function(x) sum(x > 0))
species_se <- species_sd / sqrt(species_n)
species_df <- data.frame(
  Species = colnames(species_data),
  TotalAbundance = species_totals,
  MeanAbundance = species_means,
  SE = species_se
) %>%
  arrange(desc(TotalAbundance)) %>%
  mutate(Species = factor(Species, levels = Species)) %>%
  mutate(MeanAbundance = ifelse(MeanAbundance <= 0, 0.01, MeanAbundance)) %>%
  mutate(
    FlippedMeanAbundance = MeanAbundance,
    FlippedSEMin = MeanAbundance - SE,
    FlippedSEMax = MeanAbundance + SE
  )
total_abundance_plot <- ggplot(species_df, aes(x = Species, y = TotalAbundance, fill = Species)) +
  geom_bar(stat = "identity") +
  scale_y_log10() +
  scale_fill_manual(values = species_colors) +
  theme_bw() +
  labs(y = "Log Total Abundance", x = NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  guides(fill = "none")
relative_abundance_plot <- ggplot(species_df, aes(x = Species, y = MeanAbundance, fill = Species)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(
    ymin = pmax(FlippedSEMin, 0.01),
    ymax = FlippedSEMax
  ), width = 0.2) +
  scale_fill_manual(values = species_colors) +
  scale_y_reverse(breaks = seq(0, max(species_df$MeanAbundance, na.rm = TRUE), by = 5)) +
  theme_bw() +
  labs(y = "Mean Abundance (Non-Zero, with SE)", x = NULL) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  guides(fill = "none")
combined_plot <- total_abundance_plot / relative_abundance_plot +
  plot_layout(heights = c(1, 3))
ggsave("results/Descriptive/Abundance.all.pdf", combined_plot, width = 10, height = 10)

# --- Temporal Histograms: Species Counts Over Time (All) ---
data_summary <- DATA.all %>%
  select(collectionEnd, 7:19) %>%
  pivot_longer(
    cols = -collectionEnd,
    names_to = "species",
    values_to = "count"
  ) %>%
  mutate(
    collectionEnd = as.Date(collectionEnd, format = "%d/%m/%Y"),
    interval = as.Date(cut(collectionEnd, breaks = "14 days")),
    count = as.numeric(count)
  ) %>%
  group_by(interval, species) %>%
  summarize(
    total_count = sum(count, na.rm = TRUE),
    .groups = "drop"
  )
PLOT <- ggplot(data_summary, aes(x = interval, y = total_count, fill = species)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8) +
  scale_fill_manual(values = species_colors) +
  labs(
    x = "Two-Week Intervals",
    y = "Total Count",
    title = "Histogram of Species Counts Over Two-Week Intervals"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray80", size = 0.5)
  ) +
  facet_wrap(~species, scales = "free_y")
ggsave("results/Descriptive/TemporalPresHist.all.pdf", PLOT, width = 10, height = 5)

# --- Temporal Histograms: Species Counts Over Time (Vienna) ---
data_summary <- DATA.Vienna %>%
  select(collectionEnd, 7:19) %>%
  pivot_longer(
    cols = -collectionEnd,
    names_to = "species",
    values_to = "count"
  ) %>%
  mutate(
    collectionEnd = as.Date(collectionEnd, format = "%d/%m/%Y"),
    interval = as.Date(cut(collectionEnd, breaks = "14 days")),
    count = as.numeric(count)
  ) %>%
  group_by(interval, species) %>%
  summarize(
    total_count = sum(count, na.rm = TRUE),
    .groups = "drop"
  )
PLOT <- ggplot(data_summary, aes(x = interval, y = total_count, fill = species)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8) +
  scale_fill_manual(values = species_colors) +
  labs(
    x = "Two-Week Intervals",
    y = "Total Count",
    title = "Histogram of Species Counts Over Two-Week Intervals"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray80", size = 0.5)
  ) +
  facet_wrap(~species, scales = "free_y")
ggsave("results/Descriptive/TemporalPresHist.Vienna.pdf", PLOT, width = 10, height = 5)
