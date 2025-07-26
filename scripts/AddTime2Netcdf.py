import os
import numpy as np
from netCDF4 import Dataset
from datetime import datetime, timedelta
import argparse

# Function to create a time vector based on the file name


def create_time_vector(filename, num_layers):
    if "_1d_" in filename:
        # Start from 2024/06/01 and add daily intervals
        start_date = datetime(2024, 6, 1)
        time_vector = [start_date + timedelta(days=i)
                       for i in range(num_layers)]
    elif "_1m_" in filename:
        # Start from 2024/01/01 and add monthly intervals (12 items)
        start_date = datetime(2024, 1, 1)
        time_vector = [start_date.replace(
            month=(start_date.month + i - 1) % 12 + 1) for i in range(num_layers)]
    else:
        return None  # Return None if the filename doesn't match '_1d_' or '_1m_'

    # Convert to a numpy array of datetime objects
    return np.array(time_vector, dtype='datetime64[D]')

# Function to process each NetCDF file in the directory


def process_files_in_directory(directory, verbose):
    # Loop through all files in the directory
    for filename in os.listdir(directory):
        if filename.endswith('.nc'):  # Check if the file is a NetCDF file
            file_path = os.path.join(directory, filename)

            try:
                # Open the existing NetCDF file in read/write mode
                dataset = Dataset(file_path, 'r+')

                # Try to find the first variable that likely represents the data
                data_var_name = list(dataset.variables.keys())[
                    0]  # Assuming the first variable
                data_var = dataset.variables[data_var_name]

                # Ensure the data variable has time as the first dimension
                if len(data_var.shape) < 3:  # Handle 1D or 2D variables
                    # Use the first dimension as time if it's 1D or 2D
                    num_layers = data_var.shape[0]
                else:
                    # Assuming time is the first dimension
                    num_layers = data_var.shape[0]

                # Create the time vector based on the filename
                time_vector = create_time_vector(filename, num_layers)

                if time_vector is None:
                    print(
                        f"Skipping file {filename}: Filename does not contain '_1d_' or '_1m_'")
                    continue

                # Check if the 'time' dimension exists, if not create it
                if 'time' not in dataset.dimensions:
                    # None means unlimited dimension
                    dataset.createDimension('time', None)

                # Check if the 'time' variable exists, if not create it
                if 'time' not in dataset.variables:
                    time_var = dataset.createVariable(
                        'time', np.float64, ('time',))
                    time_var[:] = time_vector

                # Save the modified NetCDF file with '_time' appended to the filename
                new_filename = filename.replace('.nc', '_time.nc')
                new_file_path = os.path.join(directory, new_filename)
                dataset.close()

                # Rename the original file to include '_time'
                os.rename(file_path, new_file_path)

                if verbose:
                    print(f"Processed and saved: {new_filename}")

            except Exception as e:
                print(f"Error processing {filename}: {e}")

# Main function to handle command-line arguments


def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(
        description='Process NetCDF files and add time dimension')

    # Arguments
    parser.add_argument('--directory', type=str, required=True,
                        help='Directory containing NetCDF files to process')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose output')

    # Parse the arguments
    args = parser.parse_args()

    # Process the files in the given directory
    process_files_in_directory(args.directory, args.verbose)


# Run the main function
if __name__ == '__main__':
    main()
