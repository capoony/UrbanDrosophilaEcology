import os
import requests
import argparse
import xarray as xr
import rioxarray
from rasterio.enums import Resampling

# Base URL for Geosphere Austria
base_url = "https://public.hub.geosphere.at/datahub/resources/spartacus-v2-"


def download_file(url, save_path):
    """Download a file from a URL and save it locally."""
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()
        with open(save_path, 'wb') as file:
            for chunk in response.iter_content(chunk_size=8192):
                file.write(chunk)
        print(f"Downloaded: {save_path}")
    except requests.exceptions.RequestException as e:
        print(f"Failed to download {url}: {e}")


def reproject_netcdf_to_wgs84(file_path):
    """Reproject the NetCDF file to WGS84 CRS and print its extent."""
    try:
        ds = xr.open_dataset(file_path)

        # Ensure CRS metadata is defined
        if not ds.rio.crs:
            ds = ds.rio.write_crs("EPSG:3416")  # Assume CRS from the dataset

        # Reproject to WGS84
        ds = ds.rio.reproject("EPSG:4326")

        # Calculate and print extent
        extent = {
            "min_lon": ds.rio.bounds()[0],
            "min_lat": ds.rio.bounds()[1],
            "max_lon": ds.rio.bounds()[2],
            "max_lat": ds.rio.bounds()[3],
        }
        print(f"Reprojected file extent: {extent}")

        # Save the reprojected file
        reprojected_path = file_path.replace(".nc", "_reprojected.nc")
        ds.to_netcdf(reprojected_path)
        print(f"Reprojected and saved: {reprojected_path}")
        return reprojected_path
    except Exception as e:
        print(f"Failed to reproject {file_path}: {e}")
        return None


# Define parameters for download
params = {
    "yearly": {
        "suffix": "YEARLY_",
        "params": ["TM", "RR", "SA"],
    },
    "monthly": {
        "suffix": "MONTHLY_",
        "params": ["TM", "RR", "SA"],
    },
    "daily": {
        "suffix": "DAILY_",
        "params": ["TX", "TN", "RR", "SA"],
    },
}

# Argument parser setup
parser = argparse.ArgumentParser(
    description="Download SPARTACUS NetCDF files.")
parser.add_argument("--output_folder", type=str, default="spartacus_data",
                    help="Base directory to save downloaded files.")
parser.add_argument("--year", type=int, default=2024,
                    help="Year for which to download data.")
args = parser.parse_args()

# Output folder and year
base_dir = args.output_folder
year = args.year
os.makedirs(base_dir, exist_ok=True)

# Loop through each time period and parameter
for period, config in params.items():
    period_dir = os.path.join(base_dir, period)
    os.makedirs(period_dir, exist_ok=True)

    for param in config["params"]:
        filename = f"SPARTACUS2-{config['suffix']}{param}_{year}.nc"
        period_suffix = "1y-" if period == "yearly" else (
            "1m-" if period == "monthly" else "1d-")
        url = f"{base_url}{period_suffix}1km/filelisting/{param}/{filename}"
        save_path = os.path.join(period_dir, filename)

        # Download the file
        download_file(url, save_path)

        # Reproject and print extent
        reproject_netcdf_to_wgs84(save_path)

print("All downloads and processing completed.")
