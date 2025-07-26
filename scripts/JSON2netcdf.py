import json
import pandas as pd
import numpy as np
import argparse
import pyproj
from scipy.interpolate import griddata
from netCDF4 import Dataset
import os


def reproject_coordinates(lon, lat, source_epsg, target_epsg):
    """
    Reproject coordinates from one CRS to another.

    Args:
        lon (float): Longitude in source CRS.
        lat (float): Latitude in source CRS.
        source_epsg (str): Source EPSG code.
        target_epsg (str): Target EPSG code.

    Returns:
        tuple: Reprojected longitude and latitude.
    """
    transformer = pyproj.Transformer.from_crs(
        source_epsg, target_epsg, always_xy=True)
    return transformer.transform(lon, lat)


def restrict_and_interpolate(df, bbox, resolution):
    """
    Restrict data to a bounding box and interpolate to a finer grid.

    Args:
        df (pd.DataFrame): DataFrame with 'longitude', 'latitude', 'date', and 'value'.
        bbox (tuple): Bounding box as (xmin, xmax, ymin, ymax).
        resolution (float): Resolution of the grid in degrees (~100m).

    Returns:
        pd.DataFrame: Interpolated DataFrame.
    """
    xmin, xmax, ymin, ymax = bbox

    # Filter data within the bounding box
    df = df[(df['longitude'] >= xmin) & (df['longitude'] <= xmax) &
            (df['latitude'] >= ymin) & (df['latitude'] <= ymax)]

    # Create the target grid
    grid_lon = np.arange(xmin, xmax, resolution)
    grid_lat = np.arange(ymin, ymax, resolution)
    grid_x, grid_y = np.meshgrid(grid_lon, grid_lat)

    interpolated_frames = []

    for date in df['date'].unique():
        df_date = df[df['date'] == date]
        interpolated_values = griddata(
            points=df_date[['longitude', 'latitude']].values,
            values=df_date['value'],
            xi=(grid_x, grid_y),
            method='linear'
        )

        interpolated_df = pd.DataFrame({
            'longitude': grid_x.ravel(),
            'latitude': grid_y.ravel(),
            'value': interpolated_values.ravel(),
            'date': date
        }).dropna()

        interpolated_frames.append(interpolated_df)

    return pd.concat(interpolated_frames, ignore_index=True)


def save_to_netcdf(df, output_file, parameter_name):
    """
    Save DataFrame to a NetCDF file.

    Args:
        df (pd.DataFrame): DataFrame with 'longitude', 'latitude', 'date', and 'value'.
        output_file (str): Path to the output NetCDF file.
        parameter_name (str): Name of the parameter to include in metadata.
    """
    with Dataset(output_file, 'w', format='NETCDF4') as ncfile:
        # Define dimensions
        lon = ncfile.createDimension('lon', len(df['longitude'].unique()))
        lat = ncfile.createDimension('lat', len(df['latitude'].unique()))
        time = ncfile.createDimension('time', len(df['date'].unique()))

        # Define variables
        lons = ncfile.createVariable('lon', 'f4', ('lon',))
        lats = ncfile.createVariable('lat', 'f4', ('lat',))
        times = ncfile.createVariable('time', 'f4', ('time',))
        values = ncfile.createVariable('value', 'f4', ('time', 'lat', 'lon',))

        # Write data
        lons[:] = sorted(df['longitude'].unique())
        lats[:] = sorted(df['latitude'].unique())
        times[:] = sorted(df['date'].unique())

        for i, date in enumerate(sorted(df['date'].unique())):
            grid = df[df['date'] == date].pivot(
                index='latitude', columns='longitude', values='value').values
            values[i, :, :] = grid

        # Add metadata
        ncfile.description = f"Interpolated data for parameter: {parameter_name}"
        ncfile.source = "Converted from GeoJSON to NetCDF"

    print(f"Saved NetCDF file to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert GeoJSON to NetCDF with reprojection and interpolation.")
    parser.add_argument('--input_file', type=str, required=True,
                        help="Path to the input GeoJSON file.")
    parser.add_argument('--output_dir', type=str, required=True,
                        help="Directory to save the output NetCDF files.")
    parser.add_argument('--bbox', type=float, nargs=4, metavar=('xmin', 'xmax', 'ymin', 'ymax'),
                        required=True, help="Bounding box as xmin xmax ymin ymax.")
    parser.add_argument('--resolution', type=float, default=0.001,
                        help="Grid resolution in degrees (~100m).")
    parser.add_argument('--source_epsg', type=str,
                        default="EPSG:4258", help="Source projection EPSG code.")
    parser.add_argument('--target_epsg', type=str,
                        default="EPSG:4326", help="Target projection EPSG code.")
    args = parser.parse_args()

    # Load GeoJSON
    with open(args.input_file, 'r') as f:
        data = json.load(f)

    input_file_name = os.path.splitext(os.path.basename(args.input_file))[0]

    # Parse features
    for param, param_data in data['features'][0]['properties']['parameters'].items():
        print(f"Processing parameter: {param}")
        records = []
        for feature in data['features']:
            lon, lat = feature['geometry']['coordinates']
            reprojected_lon, reprojected_lat = reproject_coordinates(
                lon, lat, args.source_epsg, args.target_epsg)
            for i, value in enumerate(feature['properties']['parameters'][param]['data']):
                records.append({
                    'longitude': reprojected_lon,
                    'latitude': reprojected_lat,
                    'value': value,
                    'date': i
                })

        df = pd.DataFrame(records)

        # Restrict and interpolate
        interpolated_df = restrict_and_interpolate(
            df, args.bbox, args.resolution)

        # Save to NetCDF
        output_file = f"{args.output_dir}/{input_file_name}_{param}.nc"
        save_to_netcdf(interpolated_df, output_file, param)


if __name__ == "__main__":
    main()
