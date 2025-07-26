import json
import os
import argparse
from collections import defaultdict


def parse_geojson(file_path):
    """
    Parse a GeoJSON file and return its data.

    Args:
        file_path (str): Path to the GeoJSON file.

    Returns:
        dict: Parsed GeoJSON data.
    """
    with open(file_path, 'r') as f:
        data = json.load(f)
    return data


def concatenate_geojson(input_dir, output_file):
    """
    Concatenate multiple GeoJSON files into a single file while combining values for each parameter by extending the lists of timestamps and values.

    Args:
        input_dir (str): Directory containing GeoJSON files.
        output_file (str): Path to the output GeoJSON file.
    """
    all_timestamps = set()
    combined_geojson = {
        "media_type": "application/json",
        "type": "FeatureCollection",
        "version": "v1",
        "timestamps": [],
        "features": []
    }

    # Dictionary to aggregate data by coordinates and parameters
    aggregated_data = defaultdict(lambda: defaultdict(
        lambda: {"timestamps": [], "data": []}))
    parameter_metadata = {}

    for file_name in os.listdir(input_dir):
        if file_name.endswith('.geojson') and "_1d_" in file_name:
            file_path = os.path.join(input_dir, file_name)
            print(f"Processing file: {file_path}")
            geojson_data = parse_geojson(file_path)

            # Extract top-level timestamps
            top_level_timestamps = geojson_data.get("timestamps", [])
            if not top_level_timestamps:
                raise ValueError(
                    f"No top-level timestamps found in {file_path}.")

            for feature in geojson_data.get("features", []):
                lon, lat = feature["geometry"]["coordinates"]
                parameters = feature["properties"]["parameters"]

                for param, param_data in parameters.items():
                    if param not in parameter_metadata:
                        parameter_metadata[param] = {
                            "name": param_data.get("name"),
                            "unit": param_data.get("unit")
                        }
                    values = param_data.get("data", [])

                    if len(values) != len(top_level_timestamps):
                        raise ValueError(
                            f"Mismatch in {file_path} for parameter '{param}': "
                            f"{len(values)} values vs {len(top_level_timestamps)} timestamps."
                        )

                    # Extend and align data for this coordinate and parameter
                    for ts, value in zip(top_level_timestamps, values):
                        if ts not in aggregated_data[(lon, lat)][param]["timestamps"]:
                            aggregated_data[(lon, lat)
                                            ][param]["timestamps"].append(ts)
                            aggregated_data[(lon, lat)][param]["data"].append(
                                value)

                    all_timestamps.update(top_level_timestamps)

    # Sort all timestamps
    sorted_timestamps = sorted(all_timestamps)
    combined_geojson["timestamps"] = sorted_timestamps

    # Build the combined GeoJSON
    for (lon, lat), params in aggregated_data.items():
        feature = {
            "type": "Feature",
            "geometry": {
                "type": "Point",
                "coordinates": [lon, lat]
            },
            "properties": {
                "parameters": {}
            }
        }
        for param, param_data in params.items():
            # Align data with global sorted timestamps
            aligned_data = [None] * len(sorted_timestamps)
            for ts, value in zip(param_data["timestamps"], param_data["data"]):
                index = sorted_timestamps.index(ts)
                aligned_data[index] = value

            aggregated_param_data = {
                "name": parameter_metadata[param]["name"],
                "unit": parameter_metadata[param]["unit"],
                "data": aligned_data
            }

            feature["properties"]["parameters"][param] = aggregated_param_data

        combined_geojson["features"].append(feature)

    with open(output_file, 'w') as f:
        json.dump(combined_geojson, f, indent=2)

    print(f"Concatenated GeoJSON saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Concatenate multiple GeoJSON files into a single file.")
    parser.add_argument('--input_dir', type=str, required=True,
                        help="Directory containing GeoJSON files.")
    parser.add_argument('--output_file', type=str, required=True,
                        help="Path to the output GeoJSON file.")
    args = parser.parse_args()

    concatenate_geojson(args.input_dir, args.output_file)


if __name__ == "__main__":
    main()
