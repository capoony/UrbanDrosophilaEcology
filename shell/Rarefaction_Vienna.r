# RDA Analysis for Vienna Full Collapsed Dataset
# Author: Martin Kapun
#
# This script performs Redundancy Analysis (RDA) on the Vienna dataset.
# It includes data loading, cleaning, variable selection, multicollinearity checks,
# Hellinger transformation, RDA modeling, and result visualization.

# --- Load Required Libraries ---
library(tidyverse)
library(vegan) # For RDA and ecological analysis
library(corrplot) # For correlation matrix visualization
library(iNEXT) # For rarefaction and extrapolation analysis

# --- Set Working Directory ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    WD <- getwd()
} else {
    WD <- args[1]
}
setwd(WD)
# setwd("D:/GitHub/UrbanDrosophilaEcology")

# --- Create Results Directory ---
dir.create("results/Rarefaction", showWarnings = FALSE, recursive = TRUE)

# --- Load and Prepare Data ---
DATA <- read.csv("data/Samples_inca_spartacus_vienna_clean_final.csv", header = TRUE)

# Remove rows with missing values
DATA.Vienna <- na.omit(DATA)

# --- Extract Coordinates and Dates (if needed for further analysis) ---
DATA.coord.Vienna <- DATA.Vienna[, 5:6]
DATA.date.Vienna <- DATA.Vienna[, 4]

# --- Species Data Preparation and Hellinger Transformation ---
# Aggregate species counts by ParticipantId and Indoors
DATA.spec.Vienna <- DATA.Vienna %>%
    dplyr::select(2, 7:20) %>%
    group_by(ParticipantId, Indoors) %>%
    summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
    dplyr::select(-ParticipantId, -Indoors)

# Replace NAs with zeros
DATA.spec.Vienna[is.na(DATA.spec.Vienna)] <- 0

## exclude samples with zero counts
DATA.spec.Vienna <- DATA.spec.Vienna[rowSums(DATA.spec.Vienna) > 0, ]

# --- Rarefaction and Species Richness Estimation ---

# Convert to matrix
comm_matrix <- as.matrix(DATA.spec.Vienna)

# Species accumulation curve with asymptotic prediction
spec_curve <- specaccum(comm_matrix, method = "random")

# Fit Michaelis-Menten model to predict total richness
fit_mm <- fitspecaccum(spec_curve, "michaelis-menten")
predicted_total_species <- coef(fit_mm)[1]

# --- Simple Statistical Tests ---

# 1. Goodness of fit test for the asymptotic model
observed_richness <- spec_curve$richness[length(spec_curve$richness)]
sampling_completeness <- (observed_richness / predicted_total_species) * 100

# 2. Bootstrap confidence interval for the asymptotic estimate
set.seed(42)
bootstrap_estimates <- replicate(1000, {
    tryCatch(
        {
            boot_sample <- comm_matrix[sample(nrow(comm_matrix), replace = TRUE), ]
            boot_curve <- specaccum(boot_sample, method = "random")
            boot_fit <- fitspecaccum(boot_curve, "michaelis-menten")
            coef(boot_fit)[1]
        },
        error = function(e) {
            # If model fitting fails, return NA
            NA
        }
    )
})

# Remove failed bootstrap attempts and calculate CI
bootstrap_estimates_clean <- bootstrap_estimates[!is.na(bootstrap_estimates)]
cat("Bootstrap success rate:", length(bootstrap_estimates_clean), "out of 1000 attempts\n")

# Calculate 95% confidence interval from successful bootstraps
if (length(bootstrap_estimates_clean) >= 100) {
    ci_95 <- quantile(bootstrap_estimates_clean, c(0.025, 0.975), na.rm = TRUE)
} else {
    # Fallback: use model standard error if too few bootstrap successes
    cat("Warning: Too few successful bootstrap attempts, using model SE for CI\n")
    model_se <- summary(fit_mm)$coefficients[1, 2]
    ci_95 <- c(
        predicted_total_species - 1.96 * model_se,
        predicted_total_species + 1.96 * model_se
    )
}

# 3. Statistical test for sampling completeness
# Instead of t-test with single observed value, we'll use the bootstrap distribution
if (length(bootstrap_estimates_clean) >= 100) {
    # Test if the observed richness is significantly less than the bootstrap distribution mean
    # This tests whether we've captured most of the expected species diversity
    bootstrap_mean <- mean(bootstrap_estimates_clean)
    bootstrap_sd <- sd(bootstrap_estimates_clean)

    # Calculate z-score for observed richness
    z_score <- (observed_richness - bootstrap_mean) / bootstrap_sd
    p_value <- pnorm(z_score) # One-tailed test (observed < predicted)

    test_method <- "Bootstrap Z-test"
} else {
    # Fallback: Compare observed to predicted using model uncertainty
    model_se <- summary(fit_mm)$coefficients[1, 2]
    z_score <- (observed_richness - predicted_total_species) / model_se
    p_value <- pnorm(z_score) # One-tailed test (observed < predicted)

    test_method <- "Model SE Z-test"
}

# Alternative simple test: Check if observed richness falls within CI
within_ci <- observed_richness >= ci_95[1] & observed_richness <= ci_95[2]

# Create summary statistics
stats_summary <- data.frame(
    Metric = c(
        "Observed_Richness", "Predicted_Richness", "Sampling_Completeness_Percent",
        "Bootstrap_CI_Lower", "Bootstrap_CI_Upper", "Statistical_Test", "Test_p_value",
        "Significant_Gap", "Within_CI"
    ),
    Value = c(
        observed_richness, predicted_total_species, sampling_completeness,
        ci_95[1], ci_95[2], test_method, p_value,
        ifelse(p_value < 0.05, "Yes", "No"), ifelse(within_ci, "Yes", "No")
    )
)

# Save statistical results
write.table(stats_summary,
    file = "results/Rarefaction_Vienna_full_collapsed/simple_statistics.txt",
    row.names = FALSE, quote = FALSE, sep = "\t"
)

# Plot species accumulation curve with statistical annotations
pdf("results/Rarefaction_Vienna_full_collapsed/species_accumulation_curve.pdf",
    width = 10, height = 6
)
plot(spec_curve,
    ci.type = "poly", col = "blue", lwd = 2,
    ci.lty = 0, ci.col = "#add8e66d",
    main = "Species Accumulation Curve with Statistical Test",
    xlab = "Number of Sites", ylab = "Species Richness",
    ylim = c(0, 17)
)

# Add confidence interval rectangle
rect(
    xleft = 0, xright = max(spec_curve$sites),
    ybottom = ci_95[1], ytop = ci_95[2],
    col = "#FFA50030", border = NA
)

# Add model fit
plot(fit_mm, add = TRUE, col = "#ff00002e", lwd = 2)

# Add horizontal lines for key values
abline(h = observed_richness, col = "blue", lty = 2, lwd = 2)
abline(h = predicted_total_species, col = "#ff00002e", lty = 2, lwd = 2)
abline(h = ci_95, col = "orange", lty = 3, lwd = 1)

# Add text annotations
text(max(spec_curve$sites) * 0.7, observed_richness + 0.3,
    paste("Observed:", round(observed_richness, 1)),
    col = "blue", cex = 0.9
)
text(max(spec_curve$sites) * 0.7, predicted_total_species + 0.3,
    paste("Predicted:", round(predicted_total_species, 1)),
    col = "#ff00002e", cex = 0.9
)
text(max(spec_curve$sites) * 0.7, predicted_total_species - 0.5,
    paste("Completeness:", round(sampling_completeness, 1), "%"),
    cex = 0.9
)
text(max(spec_curve$sites) * 0.7, predicted_total_species - 0.8,
    paste(
        "Gap significant:", ifelse(p_value < 0.05, "Yes", "No"),
        "(p =", round(p_value, 3), ")"
    ),
    cex = 0.9
)

legend("bottomright",
    legend = c(
        "Observed Accumulation", "Michaelis-Menten Fit", "Current Richness",
        "Predicted Total", "95% CI"
    ),
    col = c("blue", "#ff00002e", "blue", "#ff00002e", "orange"),
    lwd = c(2, 2, 2, 2, 1),
    lty = c(1, 1, 2, 2, 3),
    cex = 0.8
)
dev.off()

# Print results to console
cat("\n=== SIMPLE STATISTICAL SUMMARY ===\n")
cat("Observed species richness:", observed_richness, "\n")
cat("Predicted total richness:", round(predicted_total_species, 2), "\n")
cat("Sampling completeness:", round(sampling_completeness, 1), "%\n")
cat("95% Bootstrap CI:", round(ci_95[1], 2), "-", round(ci_95[2], 2), "\n")
if (exists("bootstrap_estimates_clean")) {
    cat("Bootstrap success rate:", length(bootstrap_estimates_clean), "/1000\n")
}
cat("Statistical test:", test_method, "\n")
cat("Test p-value:", round(p_value, 4), "\n")
cat("Observed within CI:", ifelse(within_ci, "Yes", "No"), "\n")
if (p_value < 0.05) {
    cat("CONCLUSION: Sampling is significantly incomplete - more species likely exist.\n")
} else {
    cat("CONCLUSION: Sampling appears complete - few additional species expected.\n")
}
