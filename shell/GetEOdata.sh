VmaxLat=48.5137
VminLat=47.9585
VminLon=15.9873
VmaxLon=16.8123

## set working directory
WD=/media/inter/mkapun/projects/UrbanDrosophilaEcology

### (1) INCA dataset
Rscript ${WD}/shell/GetINCAData.R

### Merge NETCDF with Samples data

python ${WD}/scripts/MergeSamplesNetCDF.py \
  --input_csv ${WD}/data/Samples.csv \
  --input_dir ${WD}/data/outputs_full \
  --output_csv ${WD}/data/Samples_inca.csv

### (2) Spartacus dataset
Rscript ${WD}/shell/GetSpartacusData.R

### (3) ViennaCube data
Rscript ${WD}/shell/GetViennaDataCubeData.r

Rscript -e '''
# Read the CSV file

setwd("${WD}")
data <- read.csv(data/Samples_inca_spartacus_vienna.csv", stringsAsFactors = FALSE)

# Rename columns based on conditions
colnames(data) <- sapply(colnames(data), function(col) {
  if (grepl("real_land_use2020", col)) {
    # Keep the substring after the last "_"
    col <- sub(".*_", "", col)
  } else if (grepl("reprojected_r\\d+_", col)) {
    # Keep the substring after "reprojected_r*_" and before "_100m"
    col <- sub(".*reprojected_r\\d+_(.*?)100m.*", "\\1", col)
  }
  
  # Remove suffix starting with "..." and trailing underscores
  col <- sub("\\.\\.\\..*", "", col)
  col <- sub("_+$", "", col)
  return(col)
})

# Write the updated CSV file
write.csv(data, "data/Samples_inca_spartacus_vienna_clean.csv", row.names = FALSE, quote = FALSE)

'''
