import os
import argparse
import pandas as pd
import numpy as np
import xarray as xr
from datetime import datetime, timedelta


def find_closest_grid(lat, lon, grid_lats, grid_lons):
    """
    Find the closest grid cell for given latitude and longitude.
    Handles 1D latitude and longitude arrays by meshing them into a 2D grid.
    """
    if grid_lats.ndim == 1 and grid_lons.ndim == 1:
        lon_grid, lat_grid = np.meshgrid(grid_lons, grid_lats)
    else:
        lat_grid, lon_grid = grid_lats, grid_lons

    distances = np.sqrt((lat_grid - lat)**2 + (lon_grid - lon)**2)
    return np.unravel_index(np.argmin(distances), distances.shape)


def process_netcdf(file_path, samples, data_type, prefix, param):
    """
    Process a NetCDF file and append the values for the closest grid cell and time index to the samples DataFrame.
    """
    ds = xr.open_dataset(file_path)

    # Dynamically handle latitude and longitude variable names
    lat_var = "latitude" if "latitude" in ds else "lat"
    lon_var = "longitude" if "longitude" in ds else "lon"
    time = "z" if "z" in ds else "time"

    grid_lats, grid_lons = np.array(ds[lat_var]), np.array(ds[lon_var])

    # Determine the start index of the time dimension
    time_start = int(ds[time][0].values) if ds[time].size > 0 else 0

    # Set the start date for daily datasets
    if data_type == "daily":
        start_date = datetime(2024, 6, 1)

    for idx, sample in samples.iterrows():
        lat, lon = sample['Latitude'], sample['Longitude']
        collection_date = datetime.strptime(
            sample['collectionEnd'], "%d/%m/%Y")

        # Find the closest grid cell
        closest_y, closest_x = find_closest_grid(
            lat, lon, grid_lats, grid_lons)

        if data_type == "daily":
            # Calculate the range of layer indices for the 14 days prior to the collection date
            days_since_start = (collection_date - start_date).days + time_start
            start_layer = max(time_start, days_since_start - 14)
            end_layer = days_since_start + 1  # Include the collection day

            # Ensure the range is within the bounds of the dataset's layers
            if time_start <= days_since_start < time_start + ds.sizes[time]:
                values = ds.isel(
                    {lat_var: closest_y, lon_var: closest_x,
                        time: slice(start_layer, end_layer)}
                ).to_array().values
                # Check if the result has more than one nested array
                if values.ndim > 1 and values.shape[0] > 1:
                    # Take only the second nested array
                    avg_value = np.nanmean(values[1])
                else:
                    avg_value = np.nanmean(values)
                samples.loc[idx, f"{prefix}_{param}_daily"] = avg_value
            else:
                samples.loc[idx, f"{prefix}_{param}_daily"] = np.nan

        elif data_type == "monthly":
            # Find the closest month index (starting from January 2024)
            month_index = (collection_date.year - 2024) * 12 + \
                collection_date.month - 1 + time_start
            if time_start <= month_index < time_start + ds.sizes[time]:
                value = ds.isel(
                    {lat_var: closest_y, lon_var: closest_x, time: month_index}
                ).to_array().values
                # Check if the result contains more than one value
                if value.size > 1:
                    value = value[1]  # Take only the second item
                samples.loc[idx, f"{prefix}_{param}_monthly"] = value.item(
                ) if value.size == 1 else np.nan
            else:
                samples.loc[idx, f"{prefix}_{param}_monthly"] = np.nan

        elif data_type == "yearly":
            # Use the only layer available in yearly datasets
            value = ds.isel(
                {lat_var: closest_y, lon_var: closest_x}
            ).mean(dim=time).to_array().values
            if value.size > 1:
                value = value[1]  # Take only the second item
            samples.loc[idx, f"{prefix}_{param}_yearly"] = value.item(
            ) if value.size == 1 else np.nan

    ds.close()
    return samples


def main():
    parser = argparse.ArgumentParser(
        description="Process NetCDF files and append values to a CSV.")
    parser.add_argument('--input_csv', type=str, required=True,
                        help="Path to the input CSV file.")
    parser.add_argument('--input_dir', type=str, required=True,
                        help="Directory containing NetCDF files.")
    parser.add_argument('--output_csv', type=str, required=True,
                        help="Path to the output CSV file.")
    args = parser.parse_args()

    samples = pd.read_csv(args.input_csv)
    netcdf_files = [f for f in os.listdir(args.input_dir) if f.endswith(".nc")]

    param_map = {
        "GL": "GlobalRadiation",
        "RH2M": "RelativeHumidity",
        "RR": "RainfallRate",
        "T2M": "Temperature2m",
        "TD2M": "DewPointTemperature",
        "UU": "WindSpeedEast",
        "VV": "WindSpeedNorth",
        "SA": "Sunshine",
        "TM": "MeanTemperature",
        "TN": "MinTemperature",
        "TX": "MaxTemperature",
        "P0": "AirPressure"
    }

    for file_name in netcdf_files:
        file_path = os.path.join(args.input_dir, file_name)
        prefix = "INCAL" if "INCAL" in file_name else "SPART"

        for param, param_desc in param_map.items():
            if param in file_name:
                if "_1d_" in file_name:
                    samples = process_netcdf(
                        file_path, samples, "daily", prefix, param_desc)
                elif "_1m_" in file_name:
                    samples = process_netcdf(
                        file_path, samples, "monthly", prefix, param_desc)
                elif "_1y_" in file_name:
                    samples = process_netcdf(
                        file_path, samples, "yearly", prefix, param_desc)

    samples.to_csv(args.output_csv, index=False)
    print(f"Processed data saved to {args.output_csv}")


if __name__ == "__main__":
    main()
