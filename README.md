# Urban Drosophila Ecology Project

## Overview

This repository contains the complete data analysis pipeline for the **Vienna City Fly Project**, a comprehensive study investigating the biodiversity and ecological patterns of Drosophila species in urban environments across Vienna, Austria. The project combines species abundance data with high-resolution environmental and climate data to understand how urban landscapes influence Drosophila community structure and distribution.

## Project Structure

```
UrbanDrosophilaEcology/
├── data/                           # Raw and processed datasets
│   ├── Samples_*.csv              # Species abundance and sampling data including Earth observation metadata
│   └── VCF_samples.xlsx           # Original input data
├── results/                       # Analysis outputs
│   ├── BioDiv_Vienna/             # Biodiversity analysis results
│   ├── RDA_all/                   # Redundancy analysis outputs
│   ├── RDA_Vienna_full/           # Vienna-specific RDA results
│   ├── Rarefaction/               # Species accumulation analysis
│   ├── Descriptive/               # Descriptive statistics
│   ├── Spatial_Autocorrelation/# Spatial autocorrelation analysis
│   ├── Temporal_Analysis_Repeated_Sites/ # Temporal trends at frequently sampled sites
│   └── SDM/                       # Species distribution models
├── scripts/                       # Python utilities
│   ├── getSpartacus.py            # SPARTACUS data retrieval
│   ├── ConcatenateDailyJSON.py    # Data concatenation
│   ├── MergeSamplesNetCDF.py      # Sample data merging
│   └── AddTime2Netcdf.py          # Temporal data processing
├── shell/                         # Main analysis scripts
│   ├── main.sh                    # Master pipeline script
│   ├── GetEOdata.sh               # Earth observation data acquisition
│   ├── BioDiv_Vienna.r            # Biodiversity analysis
│   ├── RDA_Vienna_full.r          # Redundancy analysis
│   ├── RDA_Vienna_full_collapsed.r # RDA with temporally collapsed data
│   ├── Rarefaction_Vienna.r       # Species accumulation and completeness analysis
│   ├── Descriptive.r              # Descriptive statistics
│   ├── Spatial_Autocorrelation_Analysis.r # Spatial autocorrelation and LME models
│   ├── Temporal_Analysis_Repeated_Sites.r # Temporal trends at repeated sites
│   └── SDM_Vienna.sh              # Species distribution modeling
└── README.md                      # This file
```

## Data Sources

### Species Data

- **Drosophila abundance data**: 13 species collected across Vienna sampling sites
- **Temporal coverage**: Multiple sampling dates with seasonal variation
- **Spatial coverage**: Urban gradient from city center to periphery
- **Collection metadata**: Sampling dates, coordinates, collector information

### Environmental Data

- **Vienna DataCube**: High-resolution climate and environmental layers
- **SPARTACUS data**: Urban morphology and microclimate parameters
- **INCAL data**: Meteorological variables (temperature, wind, radiation, precipitation)
- **Land use classification**: 32 categories of urban land use types
- **Temporal resolution**: Daily to monthly aggregations

## Analysis Pipeline

### 1. Data Acquisition and Preprocessing

```bash
# Master pipeline execution
bash shell/main.sh
```

**Key Steps:**

- Earth observation data retrieval (`GetEOdata.sh`)
- Climate data reprojection and alignment
- Species data Hellinger transformation
- Environmental variable standardization
- Multicollinearity assessment and variable selection

### 2. Descriptive Analysis

**Script:** `shell/Descriptive.r`

**Outputs:**

- Species abundance distributions
- Temporal patterns across seasons
- Spatial distribution maps
- Summary statistics by sampling location

### 3. Biodiversity Analysis

**Script:** `shell/BioDiv_Vienna.r`

**Methods:**

- Shannon diversity index calculation
- Simpson diversity and inverse Simpson indices
- Species richness and evenness metrics
- Principal Component Analysis (PCA) of environmental variables
- Mixed-effects models controlling for temporal and collector effects
- Non-metric multidimensional scaling (NMDS) with Bray-Curtis dissimilarity

**Outputs:**

- Diversity indices by sampling site
- PCA biplots and scree plots
- NMDS ordination plots
- Statistical significance tests

### 4. Rarefaction Analysis

**Script:** `shell/Rarefaction_Vienna.r`

**Methods:**

- Species accumulation curve analysis with asymptotic prediction
- Michaelis-Menten model fitting for total richness estimation
- Bootstrap confidence intervals (1,000 replicates) for asymptotic estimates
- Sampling completeness assessment using statistical tests
- Random sampling method for species accumulation curves

**Statistical Tests:**

- Bootstrap Z-test for sampling completeness
- 95% confidence intervals for predicted total species richness
- Goodness-of-fit evaluation for asymptotic models
- Assessment of significant gaps between observed and predicted richness

**Outputs:**

- Species accumulation curves with confidence intervals and model fits
- Statistical summary tables with completeness metrics
- Bootstrap-based uncertainty quantification
- Visual annotations showing observed vs. predicted richness patterns

### 5. Redundancy Analysis (RDA)

**Scripts:**

- `shell/RDA_Vienna_full.r` - Complete Vienna dataset
- `shell/RDA_Vienna_full_collapsed.r` - Temporally collapsed data
- `shell/RDA_AllSamples.r` - Full dataset including non-Vienna samples

**Methods:**

- Constrained ordination analysis
- Forward model selection with `ordiR2step()`
- Participant ID as conditioning variable
- Permutation tests (1,000-99,999 permutations)
- Adjusted R² calculation for explained variance

**Key Features:**

- Multicollinearity assessment with correlation matrices
- Systematic removal of redundant variables
- Environmental variable standardization
- Hellinger transformation of species data

**Outputs:**

- RDA ordination plots
- Statistical significance tests
- Variance partitioning results
- Model selection statistics

### 6. Redundancy Analysis with Collapsed Dates

**Script:** `shell/RDA_Vienna_full_collapsed.r`

**Purpose:** RDA analysis with temporally aggregated data to remove temporal variation and focus on spatial patterns.

### 7. Spatial Autocorrelation Analysis

**Script:** `shell/Spatial_Autocorrelation_Analysis.r`

**Methods:**

- Moran's I test for spatial autocorrelation in diversity indices
- Spatial correlograms across distance classes
- Linear mixed models (LME) with spatial correlation structures
- Comparison of correlation structures: Exponential, Gaussian, Spherical, Rational Quadratic
- AIC-based model selection for optimal spatial structure

**Key Features:**

- Distance-based neighbor analysis using k-nearest neighbors
- Spatial weights matrices with row standardization
- Monte-Carlo permutation tests (999 permutations) for significance
- Spatial correlation modeling with `nlme::corExp`, `corGaus`, `corSpher`, `corRatio`
- Random effects for repeated measures by location (ParticipantId)

**Outputs:**

- Moran's I statistics with significance tests for each diversity index
- Spatial correlograms showing autocorrelation patterns across distances
- LME model comparisons with AIC values and p-values in table format
- Best-fitting spatial correlation structure for each diversity metric
- Diagnostic plots for spatial residual patterns

### 8. Temporal Analysis for Repeatedly Sampled Sites

**Script:** `shell/Temporal_Analysis_Repeated_Sites.r`

**Methods:**

- Restriction to locations sampled ≥3 times (high-frequency sites)
- PERMANOVA tests for temporal patterns by month and week
- Linear mixed models with random effects for location (ParticipantId)
- Model comparison: Linear vs. Quadratic temporal trends
- Likelihood ratio tests for model selection

**Diversity Metrics Analyzed:**

- Shannon diversity index
- Species richness
- Simpson diversity
- Inverse Simpson diversity
- Pielou's evenness

**Statistical Approach:**

- Mixed-effects models: `lmer(Diversity ~ Month + (1|ParticipantId))`
- Polynomial trends: `lmer(Diversity ~ poly(Month, 2) + (1|ParticipantId))`
- AIC-based model comparison
- Chi-square tests for significance of quadratic terms

**Outputs:**

- Summary tables with AIC values for all models
- P-value tables from likelihood ratio tests
- Best model selection for each diversity index
- Temporal trend plots with population-level predictions
- PERMANOVA results for community composition changes over time

### 9. Species Distribution Modeling (SDM)

**Scripts:**

- `shell/SDM_Vienna.sh` - Random Forest approach

**Methods:**

- **Single Algorithm**: Random Forest with 500 trees
- 80/20 train-test split with stratified sampling
- 5-fold cross-validation
- Performance metrics: R², RMSE, MAE

**Outputs:**

- Habitat suitability maps (GeoTIFF format)
- Model performance statistics
- Visualization with Stadia basemaps
- Compound figures for major vs. minor species
