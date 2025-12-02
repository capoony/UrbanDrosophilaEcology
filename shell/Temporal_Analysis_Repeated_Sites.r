# Temporal_Analysis_Repeated_Sites.r
# Temporal analysis restricted to locations sampled >= 3 times
# Addressing: Temporal trends at frequently sampled locations
# Author: Martin Kapun
# Date: November 2025

cat("\n")
cat("===============================================================\n")
cat("  TEMPORAL ANALYSIS - REPEATED SITES (n >= 3)\n")
cat("===============================================================\n\n")

# =============================================================================
# PACKAGE INSTALLATION AND LOADING
# =============================================================================

cat("Checking and installing required packages...\n")

required_packages <- c("tidyverse", "vegan", "ggplot2", "gridExtra", "lubridate", "lme4", "dplyr", "mgcv")

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
library(ggplot2)
library(gridExtra)
library(lubridate)
library(lme4)
library(dplyr)

cat("All libraries loaded successfully!\n\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Pairwise PERMANOVA function
pairwise.adonis2 <- function(x, data, group_var, perm = 999) {
    cat("Running pairwise PERMANOVA...\n")
    groups <- unique(data[[group_var]])
    n_groups <- length(groups)

    if (n_groups < 2) {
        cat("Warning: Less than 2 groups found. Skipping pairwise comparisons.\n")
        return(NULL)
    }

    results <- data.frame(
        Comparison = character(),
        F_stat = numeric(),
        R2 = numeric(),
        p_value = numeric(),
        stringsAsFactors = FALSE
    )

    for (i in 1:(n_groups - 1)) {
        for (j in (i + 1):n_groups) {
            group1 <- groups[i]
            group2 <- groups[j]

            subset_idx <- data[[group_var]] %in% c(group1, group2)
            x_subset <- x[subset_idx, ]
            data_subset <- data[subset_idx, ]

            formula <- as.formula(paste0("x_subset ~ ", group_var))
            perm_result <- adonis2(formula, data = data_subset, permutations = perm, method = "bray")

            results <- rbind(results, data.frame(
                Comparison = paste(group1, "vs", group2),
                F_stat = perm_result$F[1],
                R2 = perm_result$R2[1],
                p_value = perm_result$`Pr(>F)`[1],
                stringsAsFactors = FALSE
            ))
        }
    }
    return(results)
}

# =============================================================================
# DATA LOADING AND PREPARATION
# =============================================================================

cat("Loading and preparing data...\n")

# Set working directory
WD <- "/media/inter/mkapun/projects/UrbanDrosophilaEcology"
setwd(WD)

# Create output directory
dir.create("results/Temporal_Analysis_Repeated_Sites", recursive = TRUE, showWarnings = FALSE)

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

cat("Data preparation complete.\n\n")

cat("\n=== TEMPORAL ANALYSIS - REPEATED SITES (>= 3 samples per location) ===\n")
cat("Analyzing temporal patterns at frequently sampled locations...\n\n")

# =============================================================================
# 1. IDENTIFY REPEATED SITES
# =============================================================================

# Count sampling events per location
location_counts <- DATA.Vienna %>%
    group_by(ParticipantId) %>%
    summarise(n_samples = n(), .groups = "drop") %>%
    arrange(desc(n_samples))

cat("Total unique locations:", nrow(location_counts), "\n")
cat("Locations sampled >= 3 times:", sum(location_counts$n_samples >= 3), "\n")
cat("Locations sampled >= 5 times:", sum(location_counts$n_samples >= 5), "\n")
cat("Locations sampled >= 10 times:", sum(location_counts$n_samples >= 10), "\n\n")

# Filter to locations with >= 5 samples
repeated_locations <- location_counts$ParticipantId[location_counts$n_samples >= 5]

cat("Filtering to", length(repeated_locations), "locations with >= 3 samples...\n")
cat("Sample size before filtering:", nrow(DATA.Vienna), "\n")

# Subset data
DATA.repeated_idx <- DATA.Vienna$ParticipantId %in% repeated_locations
DATA.Vienna.repeated <- DATA.Vienna[DATA.repeated_idx, ]
DATA.spec.repeated <- DATA.spec[DATA.repeated_idx, ]

# Recalculate diversity indices for repeated sites
shannon_div_repeated <- diversity(DATA.spec.repeated, index = "shannon")
simpson_div_repeated <- diversity(DATA.spec.repeated, index = "simpson")
invsimpson_div_repeated <- diversity(DATA.spec.repeated, index = "invsimpson")
richness_repeated <- specnumber(DATA.spec.repeated)
evenness_repeated <- shannon_div_repeated / log(richness_repeated)

cat("Sample size after filtering:", nrow(DATA.Vienna.repeated), "\n")
cat("Samples retained:", round(100 * nrow(DATA.Vienna.repeated) / nrow(DATA.Vienna), 1), "%\n\n")

# =============================================================================
# 2. PERMANOVA ANALYSES - REPEATED SITES
# =============================================================================

sink("results/Temporal_Analysis_Repeated_Sites/temporal_analysis_repeated_sites.txt")

cat("=== TEMPORAL ANALYSIS - REPEATED SITES (>= 3 samples per location) ===\n\n")
cat("Total locations:", length(repeated_locations), "\n")
cat("Total samples:", nrow(DATA.Vienna.repeated), "\n")
cat("Proportion of full dataset:", round(100 * nrow(DATA.Vienna.repeated) / nrow(DATA.Vienna), 1), "%\n\n")

# Summary of sampling frequency
cat("--- Sampling Frequency Distribution ---\n")
freq_table <- table(location_counts$n_samples[location_counts$n_samples >= 3])
print(freq_table)
cat("\n")

# 2.1 PERMANOVA by Month (repeated sites)
cat("\n--- PERMANOVA: Community Composition by Month (Repeated Sites) ---\n")
month_permanova_repeated <- adonis2(DATA.spec.repeated ~ Month,
    data = DATA.Vienna.repeated,
    permutations = 9999, method = "bray"
)
print(month_permanova_repeated)

# 2.2 PERMANOVA by Week (repeated sites)
cat("\n--- PERMANOVA: Community Composition by Week (Repeated Sites) ---\n")
week_permanova_repeated <- adonis2(DATA.spec.repeated ~ Week,
    data = DATA.Vienna.repeated,
    permutations = 9999, method = "bray"
)
print(week_permanova_repeated)

# 2.3 PERMANOVA with Month + Location (accounting for spatial pseudoreplication)
cat("\n--- PERMANOVA: Month + Location (Repeated Sites) ---\n")
month_location_permanova <- adonis2(DATA.spec.repeated ~ Month + ParticipantId,
    data = DATA.Vienna.repeated,
    permutations = 9999, method = "bray"
)
print(month_location_permanova)

# 2.4 PERMANOVA with Month + Indoors interaction
cat("\n--- PERMANOVA: Month * Indoors (Repeated Sites) ---\n")
month_indoor_permanova_repeated <- adonis2(DATA.spec.repeated ~ Month * Indoors,
    data = DATA.Vienna.repeated,
    permutations = 9999, method = "bray"
)
print(month_indoor_permanova_repeated)

# 2.5 Pairwise PERMANOVA between months (repeated sites)
cat("\n--- Pairwise PERMANOVA between Months (Repeated Sites) ---\n")
if (exists("pairwise.adonis2")) {
    pairwise_month_repeated <- pairwise.adonis2(DATA.spec.repeated,
        data = DATA.Vienna.repeated,
        group_var = "MonthName"
    )
    print(pairwise_month_repeated)
} else {
    cat("Warning: pairwise.adonis2 function not found. Skipping pairwise comparisons.\n")
}

sink()

# =============================================================================
# 3. NMDS ORDINATION - REPEATED SITES
# =============================================================================

cat("Running NMDS ordination for repeated sites...\n")
nmds_repeated <- metaMDS(DATA.spec.repeated,
    distance = "bray", k = 2, trymax = 200,
    autotransform = TRUE, trace = FALSE
)

# Check convergence
if (nmds_repeated$converged) {
    cat("NMDS converged successfully. Stress =", round(nmds_repeated$stress, 3), "\n")
} else {
    cat("Warning: NMDS did not fully converge. Stress =", round(nmds_repeated$stress, 3), "\n")
}

# Extract scores and add metadata
nmds_scores_repeated <- as.data.frame(scores(nmds_repeated, display = "sites"))
nmds_scores_repeated$Month <- DATA.Vienna.repeated$Month
nmds_scores_repeated$MonthName <- DATA.Vienna.repeated$MonthName
nmds_scores_repeated$Week <- DATA.Vienna.repeated$Week
nmds_scores_repeated$Indoors <- DATA.Vienna.repeated$Indoors
nmds_scores_repeated$ParticipantId <- DATA.Vienna.repeated$ParticipantId

# Plot NMDS by Month
p_nmds_month_repeated <- ggplot(nmds_scores_repeated, aes(x = NMDS1, y = NMDS2, color = MonthName)) +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(level = 0.95, linetype = 2) +
    labs(
        title = "NMDS - Repeated Sites Only (>= 3 samples)",
        subtitle = paste0(
            "Stress = ", round(nmds_repeated$stress, 3),
            " | n = ", nrow(DATA.Vienna.repeated), " samples from ",
            length(repeated_locations), " locations"
        )
    ) +
    theme_bw() +
    scale_color_viridis_d() +
    theme(legend.position = "right")

ggsave("results/Temporal_Analysis_Repeated_Sites/NMDS_by_Month_Repeated_Sites.pdf",
    p_nmds_month_repeated,
    width = 10, height = 7
)
ggsave("results/Temporal_Analysis_Repeated_Sites/NMDS_by_Month_Repeated_Sites.png",
    p_nmds_month_repeated,
    width = 10, height = 7, dpi = 300
)

# Plot by location to show repeated sampling
p_nmds_location <- ggplot(nmds_scores_repeated, aes(x = NMDS1, y = NMDS2, color = ParticipantId)) +
    geom_point(size = 2, alpha = 0.6) +
    geom_path(aes(group = ParticipantId), alpha = 0.3) +
    labs(
        title = "NMDS - Temporal Trajectories by Location",
        subtitle = "Lines connect samples from same location over time"
    ) +
    theme_bw() +
    theme(legend.position = "none")

ggsave("results/Temporal_Analysis_Repeated_Sites/NMDS_Trajectories_by_Location.pdf",
    p_nmds_location,
    width = 10, height = 7
)
ggsave("results/Temporal_Analysis_Repeated_Sites/NMDS_Trajectories_by_Location.png",
    p_nmds_location,
    width = 10, height = 7, dpi = 300
)

cat("NMDS plots saved.\n")

# =============================================================================
# 4. TEMPORAL TRENDS IN DIVERSITY - REPEATED SITES
# =============================================================================

div_temporal_repeated <- data.frame(
    Month = DATA.Vienna.repeated$Month,
    MonthName = DATA.Vienna.repeated$MonthName,
    Week = DATA.Vienna.repeated$Week,
    Shannon = shannon_div_repeated,
    Simpson = simpson_div_repeated,
    InvSimpson = invsimpson_div_repeated,
    Richness = richness_repeated,
    Eveness = evenness_repeated,
    Indoors = DATA.Vienna.repeated$Indoors,
    ParticipantId = DATA.Vienna.repeated$ParticipantId
)

# Linear Mixed Models for temporal trends (accounting for repeated measures by location)
sink("results/Temporal_Analysis_Repeated_Sites/temporal_analysis_repeated_sites.txt", append = TRUE)
cat("\n=== TEMPORAL TRENDS IN DIVERSITY (REPEATED SITES) ===\n")
cat("\nAnalysis based on Linear Mixed Models with location as random effect.\n")
cat("Testing linear and quadratic (polynomial) temporal trends.\n\n")

if (require(lme4, quietly = TRUE)) {
    cat("=== LINEAR MIXED MODELS: Accounting for Location as Random Effect ===\n")
    
    # Fit all mixed models
    cat("\n--- Fitting Linear Mixed Models ---\n")
    lmer_shannon_linear <- lmer(Shannon ~ Month + (1 | ParticipantId), data = div_temporal_repeated)
    lmer_richness_linear <- lmer(Richness ~ Month + (1 | ParticipantId), data = div_temporal_repeated)
    lmer_simpson_linear <- lmer(Simpson ~ Month + (1 | ParticipantId), data = div_temporal_repeated)
    lmer_invsimpson_linear <- lmer(InvSimpson ~ Month + (1 | ParticipantId), data = div_temporal_repeated)
    lmer_eveness_linear <- lmer(Eveness ~ Month + (1 | ParticipantId), data = div_temporal_repeated)
    
    cat("--- Fitting Quadratic Mixed Models ---\n")
    lmer_shannon_poly <- lmer(Shannon ~ poly(Month, 2) + (1 | ParticipantId), data = div_temporal_repeated)
    lmer_richness_poly <- lmer(Richness ~ poly(Month, 2) + (1 | ParticipantId), data = div_temporal_repeated)
    lmer_simpson_poly <- lmer(Simpson ~ poly(Month, 2) + (1 | ParticipantId), data = div_temporal_repeated)
    lmer_invsimpson_poly <- lmer(InvSimpson ~ poly(Month, 2) + (1 | ParticipantId), data = div_temporal_repeated)
    lmer_eveness_poly <- lmer(Eveness ~ poly(Month, 2) + (1 | ParticipantId), data = div_temporal_repeated)

    # Detailed model summaries
    cat("\n\n--- SHANNON DIVERSITY MODELS ---\n")
    cat("\n1. Linear Model - Shannon ~ Month + (1|ParticipantId):\n")
    print(summary(lmer_shannon_linear))
    cat("AIC:", AIC(lmer_shannon_linear), "\n")

    cat("\n2. Quadratic Model - Shannon ~ poly(Month, 2) + (1|ParticipantId):\n")
    print(summary(lmer_shannon_poly))
    cat("AIC:", AIC(lmer_shannon_poly), "\n")

    cat("\n\n--- SPECIES RICHNESS MODELS ---\n")
    cat("\n3. Linear Model - Richness ~ Month + (1|ParticipantId):\n")
    print(summary(lmer_richness_linear))
    cat("AIC:", AIC(lmer_richness_linear), "\n")

    cat("\n4. Quadratic Model - Richness ~ poly(Month, 2) + (1|ParticipantId):\n")
    print(summary(lmer_richness_poly))
    cat("AIC:", AIC(lmer_richness_poly), "\n")

    cat("\n\n--- SIMPSON DIVERSITY MODELS ---\n")
    cat("\n5. Linear Model - Simpson ~ Month + (1|ParticipantId):\n")
    print(summary(lmer_simpson_linear))
    cat("AIC:", AIC(lmer_simpson_linear), "\n")

    cat("\n6. Quadratic Model - Simpson ~ poly(Month, 2) + (1|ParticipantId):\n")
    print(summary(lmer_simpson_poly))
    cat("AIC:", AIC(lmer_simpson_poly), "\n")

    cat("\n\n--- INVERSE SIMPSON DIVERSITY MODELS ---\n")
    cat("\n7. Linear Model - InvSimpson ~ Month + (1|ParticipantId):\n")
    print(summary(lmer_invsimpson_linear))
    cat("AIC:", AIC(lmer_invsimpson_linear), "\n")

    cat("\n8. Quadratic Model - InvSimpson ~ poly(Month, 2) + (1|ParticipantId):\n")
    print(summary(lmer_invsimpson_poly))
    cat("AIC:", AIC(lmer_invsimpson_poly), "\n")

    cat("\n\n--- EVENESS MODELS ---\n")
    cat("\n9. Linear Model - Eveness ~ Month + (1|ParticipantId):\n")
    print(summary(lmer_eveness_linear))
    cat("AIC:", AIC(lmer_eveness_linear), "\n")

    cat("\n10. Quadratic Model - Eveness ~ poly(Month, 2) + (1|ParticipantId):\n")
    print(summary(lmer_eveness_poly))
    cat("AIC:", AIC(lmer_eveness_poly), "\n")

    # ANOVA comparisons
    cat("\n\n=== MODEL COMPARISONS (Likelihood Ratio Tests) ===\n")
    
    cat("\nShannon: Linear vs Quadratic\n")
    anova_shannon <- anova(lmer_shannon_linear, lmer_shannon_poly)
    print(anova_shannon)

    cat("\nRichness: Linear vs Quadratic\n")
    anova_richness <- anova(lmer_richness_linear, lmer_richness_poly)
    print(anova_richness)

    cat("\nSimpson: Linear vs Quadratic\n")
    anova_simpson <- anova(lmer_simpson_linear, lmer_simpson_poly)
    print(anova_simpson)

    cat("\nInvSimpson: Linear vs Quadratic\n")
    anova_invsimpson <- anova(lmer_invsimpson_linear, lmer_invsimpson_poly)
    print(anova_invsimpson)

    cat("\nEveness: Linear vs Quadratic\n")
    anova_eveness <- anova(lmer_eveness_linear, lmer_eveness_poly)
    print(anova_eveness)

    # Create summary tables
    cat("\n\n=== SUMMARY TABLES ===\n")
    
    # AIC Table
    cat("\n--- TABLE 1: AIC Values ---\n")
    aic_table <- data.frame(
        Diversity_Index = c("Shannon", "Richness", "Simpson", "InvSimpson", "Eveness"),
        Linear_AIC = c(AIC(lmer_shannon_linear), AIC(lmer_richness_linear), 
                       AIC(lmer_simpson_linear), AIC(lmer_invsimpson_linear), 
                       AIC(lmer_eveness_linear)),
        Quadratic_AIC = c(AIC(lmer_shannon_poly), AIC(lmer_richness_poly), 
                          AIC(lmer_simpson_poly), AIC(lmer_invsimpson_poly), 
                          AIC(lmer_eveness_poly))
    )
    aic_table$Delta_AIC <- aic_table$Linear_AIC - aic_table$Quadratic_AIC
    aic_table$Best_Model <- ifelse(aic_table$Delta_AIC > 0, "Quadratic", "Linear")
    print(aic_table, row.names = FALSE)
    
    # P-value Table
    cat("\n\n--- TABLE 2: P-values (Likelihood Ratio Tests) ---\n")
    pval_table <- data.frame(
        Diversity_Index = c("Shannon", "Richness", "Simpson", "InvSimpson", "Eveness"),
        P_value = c(
            anova_shannon$`Pr(>Chisq)`[2],
            anova_richness$`Pr(>Chisq)`[2],
            anova_simpson$`Pr(>Chisq)`[2],
            anova_invsimpson$`Pr(>Chisq)`[2],
            anova_eveness$`Pr(>Chisq)`[2]
        ),
        Chi_square = c(
            anova_shannon$Chisq[2],
            anova_richness$Chisq[2],
            anova_simpson$Chisq[2],
            anova_invsimpson$Chisq[2],
            anova_eveness$Chisq[2]
        ),
        Df = rep(1, 5)
    )
    pval_table$Significance <- ifelse(pval_table$P_value < 0.001, "***", 
                                      ifelse(pval_table$P_value < 0.01, "**",
                                            ifelse(pval_table$P_value < 0.05, "*",
                                                  ifelse(pval_table$P_value < 0.1, ".", "NS"))))
    print(pval_table, row.names = FALSE)
    
    cat("\n--- TABLE 3: Best Models Summary ---\n")
    best_models <- data.frame(
        Diversity_Index = c("Shannon", "Richness", "Simpson", "InvSimpson", "Eveness"),
        Best_Model = aic_table$Best_Model,
        AIC = ifelse(aic_table$Best_Model == "Quadratic", 
                    aic_table$Quadratic_AIC, aic_table$Linear_AIC),
        P_value = pval_table$P_value,
        Significance = pval_table$Significance
    )
    print(best_models, row.names = FALSE)
    
} else {
    cat("\nERROR: lme4 package not available for mixed models.\n")
}

sink()

# =============================================================================
# 5. DIVERSITY PLOTS - REPEATED SITES
# =============================================================================

# Generate predictions from mixed models for plotting
if (require(lme4, quietly = TRUE)) {
    # Prepare prediction data for smooth curves
    month_seq <- seq(min(div_temporal_repeated$Month), max(div_temporal_repeated$Month), length.out = 100)
    pred_data <- data.frame(Month = month_seq)
    
    # For mixed models, we predict at the population level (random effects = 0)
    # Using a "typical" location for visualization
    pred_data_mixed <- data.frame(
        Month = month_seq,
        ParticipantId = div_temporal_repeated$ParticipantId[1]  # Use first location as reference
    )
    
    # Generate predictions from mixed models
    pred_data$Shannon_Linear <- predict(lmer_shannon_linear, newdata = pred_data_mixed, re.form = NA)
    pred_data$Shannon_Poly <- predict(lmer_shannon_poly, newdata = pred_data_mixed, re.form = NA)
    pred_data$Richness_Linear <- predict(lmer_richness_linear, newdata = pred_data_mixed, re.form = NA)
    pred_data$Richness_Poly <- predict(lmer_richness_poly, newdata = pred_data_mixed, re.form = NA)
    pred_data$Simpson_Linear <- predict(lmer_simpson_linear, newdata = pred_data_mixed, re.form = NA)
    pred_data$Simpson_Poly <- predict(lmer_simpson_poly, newdata = pred_data_mixed, re.form = NA)
    pred_data$InvSimpson_Linear <- predict(lmer_invsimpson_linear, newdata = pred_data_mixed, re.form = NA)
    pred_data$InvSimpson_Poly <- predict(lmer_invsimpson_poly, newdata = pred_data_mixed, re.form = NA)
    pred_data$Eveness_Linear <- predict(lmer_eveness_linear, newdata = pred_data_mixed, re.form = NA)
    pred_data$Eveness_Poly <- predict(lmer_eveness_poly, newdata = pred_data_mixed, re.form = NA)
}

# Shannon diversity by month
p_shannon_month_repeated <- ggplot(div_temporal_repeated, aes(x = Month, y = Shannon)) +
    geom_boxplot(aes(group = Month), fill = "lightblue", alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    stat_summary(fun = mean, geom = "point", color = "black", size = 3) +
    geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = TRUE, color = "red", size = 1.2, alpha = 0.2) +
    scale_x_continuous(breaks = unique(div_temporal_repeated$Month), 
                       labels = unique(div_temporal_repeated$MonthName)[order(unique(div_temporal_repeated$Month))]) +
    labs(
        title = "Shannon Diversity by Month - Repeated Sites (with Quadratic Fit)",
        subtitle = paste0("n = ", nrow(div_temporal_repeated), " samples"),
        x = "Month",
        y = "Shannon Index"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Species richness by month
p_richness_month_repeated <- ggplot(div_temporal_repeated, aes(x = Month, y = Richness)) +
    geom_boxplot(aes(group = Month), fill = "lightgreen", alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    stat_summary(fun = mean, geom = "point", color = "black", size = 3) +
    geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = TRUE, color = "red", size = 1.2, alpha = 0.2) +
    scale_x_continuous(breaks = unique(div_temporal_repeated$Month), 
                       labels = unique(div_temporal_repeated$MonthName)[order(unique(div_temporal_repeated$Month))]) +
    labs(
        title = "Species Richness by Month - Repeated Sites (with Quadratic Fit)",
        subtitle = paste0("n = ", nrow(div_temporal_repeated), " samples"),
        x = "Month",
        y = "Species Richness"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Simpson diversity by month
p_simpson_month_repeated <- ggplot(div_temporal_repeated, aes(x = Month, y = Simpson)) +
    geom_boxplot(aes(group = Month), fill = "lightyellow", alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    stat_summary(fun = mean, geom = "point", color = "black", size = 3) +
    geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = TRUE, color = "red", size = 1.2, alpha = 0.2) +
    scale_x_continuous(breaks = unique(div_temporal_repeated$Month), 
                       labels = unique(div_temporal_repeated$MonthName)[order(unique(div_temporal_repeated$Month))]) +
    labs(
        title = "Simpson Diversity by Month - Repeated Sites (with Quadratic Fit)",
        subtitle = paste0("n = ", nrow(div_temporal_repeated), " samples"),
        x = "Month",
        y = "Simpson Index"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Inverse Simpson diversity by month
p_invsimpson_month_repeated <- ggplot(div_temporal_repeated, aes(x = Month, y = InvSimpson)) +
    geom_boxplot(aes(group = Month), fill = "lightpink", alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    stat_summary(fun = mean, geom = "point", color = "black", size = 3) +
    geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = TRUE, color = "red", size = 1.2, alpha = 0.2) +
    scale_x_continuous(breaks = unique(div_temporal_repeated$Month), 
                       labels = unique(div_temporal_repeated$MonthName)[order(unique(div_temporal_repeated$Month))]) +
    labs(
        title = "Inverse Simpson Diversity by Month - Repeated Sites (with Quadratic Fit)",
        subtitle = paste0("n = ", nrow(div_temporal_repeated), " samples"),
        x = "Month",
        y = "Inverse Simpson Index"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Eveness by month
p_eveness_month_repeated <- ggplot(div_temporal_repeated, aes(x = Month, y = Eveness)) +
    geom_boxplot(aes(group = Month), fill = "lightcyan", alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    stat_summary(fun = mean, geom = "point", color = "black", size = 3) +
    geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = TRUE, color = "red", size = 1.2, alpha = 0.2) +
    scale_x_continuous(breaks = unique(div_temporal_repeated$Month), 
                       labels = unique(div_temporal_repeated$MonthName)[order(unique(div_temporal_repeated$Month))]) +
    labs(
        title = "Eveness by Month - Repeated Sites (with Quadratic Fit)",
        subtitle = paste0("n = ", nrow(div_temporal_repeated), " samples"),
        x = "Month",
        y = "Eveness"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine all boxplots
p_temporal_combined_repeated <- grid.arrange(p_shannon_month_repeated, p_richness_month_repeated,
    p_simpson_month_repeated, p_invsimpson_month_repeated,
    p_eveness_month_repeated,
    ncol = 3, nrow = 2
)

ggsave("results/Temporal_Analysis_Repeated_Sites/Diversity_by_Month_Repeated_Sites.pdf",
    p_temporal_combined_repeated,
    width = 18, height = 10
)
ggsave("results/Temporal_Analysis_Repeated_Sites/Diversity_by_Month_Repeated_Sites.png",
    p_temporal_combined_repeated,
    width = 18, height = 10, dpi = 300
)

# --- Mixed Model Fit Comparison Plots ---
if (require(lme4, quietly = TRUE)) {
    # Shannon diversity with mixed model fits
    p_shannon_fits <- ggplot(div_temporal_repeated, aes(x = Month, y = Shannon)) +
        geom_point(alpha = 0.5, size = 2) +
        geom_line(data = pred_data, aes(y = Shannon_Linear, color = "Linear"), linetype = "dashed", size = 1) +
        geom_line(data = pred_data, aes(y = Shannon_Poly, color = "Quadratic"), size = 1.2) +
        scale_color_manual(
            name = "Mixed Model",
            values = c("Linear" = "blue", "Quadratic" = "red"),
            labels = c(
                paste0("Linear (AIC=", round(AIC(lmer_shannon_linear), 1), ")"),
                paste0("Quadratic (AIC=", round(AIC(lmer_shannon_poly), 1), ")")
            )
        ) +
        labs(
            title = "Shannon Diversity: Mixed Model Comparison",
            subtitle = "Population-level predictions (random effects = 0)",
            x = "Month (1=Jan, 12=Dec)", y = "Shannon Index"
        ) +
        theme_bw() +
        theme(legend.position = "right")

    # Species richness with mixed model fits
    p_richness_fits <- ggplot(div_temporal_repeated, aes(x = Month, y = Richness)) +
        geom_point(alpha = 0.5, size = 2) +
        geom_line(data = pred_data, aes(y = Richness_Linear, color = "Linear"), linetype = "dashed", size = 1) +
        geom_line(data = pred_data, aes(y = Richness_Poly, color = "Quadratic"), size = 1.2) +
        scale_color_manual(
            name = "Mixed Model",
            values = c("Linear" = "blue", "Quadratic" = "red"),
            labels = c(
                paste0("Linear (AIC=", round(AIC(lmer_richness_linear), 1), ")"),
                paste0("Quadratic (AIC=", round(AIC(lmer_richness_poly), 1), ")")
            )
        ) +
        labs(
            title = "Species Richness: Mixed Model Comparison",
            subtitle = "Population-level predictions (random effects = 0)",
            x = "Month (1=Jan, 12=Dec)", y = "Species Richness"
        ) +
        theme_bw() +
        theme(legend.position = "right")

    # Save model comparison plots
    p_fits_combined <- grid.arrange(p_shannon_fits, p_richness_fits, ncol = 2)

    ggsave("results/Temporal_Analysis_Repeated_Sites/Temporal_Trend_Model_Comparison.pdf",
        p_fits_combined,
        width = 16, height = 6
    )
    ggsave("results/Temporal_Analysis_Repeated_Sites/Temporal_Trend_Model_Comparison.png",
        p_fits_combined,
        width = 16, height = 6, dpi = 300
    )

    # --- Additional Diversity Indices Mixed Model Comparison Plots ---

    # Simpson diversity model fits
    p_simpson_fits <- ggplot(div_temporal_repeated, aes(x = Month, y = Simpson)) +
        geom_point(alpha = 0.5, size = 2) +
        geom_line(data = pred_data, aes(y = Simpson_Linear, color = "Linear"), linetype = "dashed", size = 1) +
        geom_line(data = pred_data, aes(y = Simpson_Poly, color = "Quadratic"), size = 1.2) +
        scale_color_manual(
            name = "Mixed Model",
            values = c("Linear" = "blue", "Quadratic" = "red"),
            labels = c(
                paste0("Linear (AIC=", round(AIC(lmer_simpson_linear), 1), ")"),
                paste0("Quadratic (AIC=", round(AIC(lmer_simpson_poly), 1), ")")
            )
        ) +
        labs(
            title = "Simpson Diversity: Mixed Model Comparison",
            x = "Month (1=Jan, 12=Dec)", y = "Simpson Index"
        ) +
        theme_bw() +
        theme(legend.position = "right")

    # InvSimpson diversity model fits
    p_invsimpson_fits <- ggplot(div_temporal_repeated, aes(x = Month, y = InvSimpson)) +
        geom_point(alpha = 0.5, size = 2) +
        geom_line(data = pred_data, aes(y = InvSimpson_Linear, color = "Linear"), linetype = "dashed", size = 1) +
        geom_line(data = pred_data, aes(y = InvSimpson_Poly, color = "Quadratic"), size = 1.2) +
        scale_color_manual(
            name = "Mixed Model",
            values = c("Linear" = "blue", "Quadratic" = "red"),
            labels = c(
                paste0("Linear (AIC=", round(AIC(lmer_invsimpson_linear), 1), ")"),
                paste0("Quadratic (AIC=", round(AIC(lmer_invsimpson_poly), 1), ")")
            )
        ) +
        labs(
            title = "Inverse Simpson: Mixed Model Comparison",
            x = "Month (1=Jan, 12=Dec)", y = "Inverse Simpson Index"
        ) +
        theme_bw() +
        theme(legend.position = "right")

    # Eveness model fits
    p_eveness_fits <- ggplot(div_temporal_repeated, aes(x = Month, y = Eveness)) +
        geom_point(alpha = 0.5, size = 2) +
        geom_line(data = pred_data, aes(y = Eveness_Linear, color = "Linear"), linetype = "dashed", size = 1) +
        geom_line(data = pred_data, aes(y = Eveness_Poly, color = "Quadratic"), size = 1.2) +
        scale_color_manual(
            name = "Mixed Model",
            values = c("Linear" = "blue", "Quadratic" = "red"),
            labels = c(
                paste0("Linear (AIC=", round(AIC(lmer_eveness_linear), 1), ")"),
                paste0("Quadratic (AIC=", round(AIC(lmer_eveness_poly), 1), ")")
            )
        ) +
        labs(
            title = "Eveness: Mixed Model Comparison",
            x = "Month (1=Jan, 12=Dec)", y = "Eveness"
        ) +
        theme_bw() +
        theme(legend.position = "right")

    # Save additional model comparison plots (2x2 grid with 3 plots)
    p_additional_fits <- grid.arrange(p_simpson_fits, p_invsimpson_fits, p_eveness_fits, ncol = 2, nrow = 2)

    ggsave("results/Temporal_Analysis_Repeated_Sites/Additional_Indices_Model_Comparison.pdf",
        p_additional_fits,
        width = 16, height = 12
    )
    ggsave("results/Temporal_Analysis_Repeated_Sites/Additional_Indices_Model_Comparison.png",
        p_additional_fits,
        width = 16, height = 12, dpi = 300
    )
}

# =============================================================================
# 6. LOCATION-SPECIFIC TEMPORAL TRENDS
# =============================================================================

# Plot temporal trends for each frequently sampled location
top_locations <- location_counts$ParticipantId[location_counts$n_samples >= 5]

if (length(top_locations) > 0) {
    cat("\nPlotting temporal trends for locations with >= 5 samples (n =", length(top_locations), ")...\n")

    div_top_locations <- div_temporal_repeated %>%
        filter(ParticipantId %in% top_locations)

    # Add quadratic trend lines by location
    p_location_trends <- ggplot(div_top_locations, aes(x = Month, y = Shannon, color = ParticipantId)) +
        geom_point(alpha = 0.6, size = 2) +
        geom_line(alpha = 0.4) +
        geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = TRUE, alpha = 0.2, size = 0.8) +
        facet_wrap(~ParticipantId, scales = "free_y", ncol = 4) +
        labs(
            title = "Shannon Diversity Temporal Trends by Location (with Quadratic Fits)",
            subtitle = paste0("Locations sampled >= 5 times (n = ", length(top_locations), ")"),
            x = "Month",
            y = "Shannon Index"
        ) +
        theme_bw() +
        theme(legend.position = "none",
              strip.text = element_text(size = 8))

    ggsave("results/Temporal_Analysis_Repeated_Sites/Diversity_Trends_by_Location.pdf",
        p_location_trends,
        width = 14, height = 10
    )
    ggsave("results/Temporal_Analysis_Repeated_Sites/Diversity_Trends_by_Location.png",
        p_location_trends,
        width = 14, height = 10, dpi = 300
    )
}

cat("Diversity plots saved.\n")

# =============================================================================
# 7. SUMMARY STATISTICS
# =============================================================================

sink("results/Temporal_Analysis_Repeated_Sites/temporal_analysis_repeated_sites.txt", append = TRUE)

cat("\n\n=== SUMMARY STATISTICS ===\n\n")

cat("--- Sampling Summary ---\n")
cat("Total locations with >= 3 samples:", length(repeated_locations), "\n")
cat("Total samples from repeated sites:", nrow(DATA.Vienna.repeated), "\n")
cat("Proportion of full dataset:", round(100 * nrow(DATA.Vienna.repeated) / nrow(DATA.Vienna), 1), "%\n")
cat("Locations with >= 5 samples:", sum(location_counts$n_samples >= 5), "\n")
cat("Locations with >= 10 samples:", sum(location_counts$n_samples >= 10), "\n\n")

cat("--- Diversity Statistics (Repeated Sites) ---\n")
cat(
    "Shannon diversity: mean =", round(mean(shannon_div_repeated), 3),
    ", SD =", round(sd(shannon_div_repeated), 3), "\n"
)
cat(
    "Species richness: mean =", round(mean(richness_repeated), 2),
    ", SD =", round(sd(richness_repeated), 2), "\n\n"
)

# Calculate full dataset diversity for comparison
shannon_div_full <- diversity(DATA.spec, index = "shannon")
richness_full <- specnumber(DATA.spec)

# Compare to full dataset
cat("--- Comparison to Full Dataset ---\n")
cat(
    "Full dataset Shannon: mean =", round(mean(shannon_div_full), 3),
    ", SD =", round(sd(shannon_div_full), 3), "\n"
)
cat(
    "Repeated sites Shannon: mean =", round(mean(shannon_div_repeated), 3),
    ", SD =", round(sd(shannon_div_repeated), 3), "\n"
)
cat("Difference:", round(mean(shannon_div_repeated) - mean(shannon_div_full), 3), "\n\n")

cat(
    "Full dataset Richness: mean =", round(mean(richness_full), 2),
    ", SD =", round(sd(richness_full), 2), "\n"
)
cat(
    "Repeated sites Richness: mean =", round(mean(richness_repeated), 2),
    ", SD =", round(sd(richness_repeated), 2), "\n"
)
cat("Difference:", round(mean(richness_repeated) - mean(richness_full), 2), "\n\n")

sink()

cat("\n=== TEMPORAL ANALYSIS (REPEATED SITES) COMPLETE ===\n")
cat("Results saved to: results/Temporal_Analysis_Repeated_Sites/temporal_analysis_repeated_sites.txt\n")
cat("Plots saved to: results/Temporal_Analysis_Repeated_Sites/\n\n")
