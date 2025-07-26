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
│   ├── Descriptive/               # Descriptive statistics
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
│   ├── Descriptive.r              # Descriptive statistics
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

### 4. Redundancy Analysis (RDA)

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

### 5. Species Distribution Modeling (SDM)

**Scripts:**

- `shell/SDM_Vienna.sh` - Random Forest approach

**Methods:**

- **Single Algorithm (Standard)**: Random Forest with 500 trees
- 80/20 train-test split with stratified sampling
- 5-fold cross-validation
- Performance metrics: R², RMSE, MAE

**Outputs:**

- Habitat suitability maps (GeoTIFF format)
- Model performance statistics
- Visualization with Stadia basemaps
- Compound figures for major vs. minor species
