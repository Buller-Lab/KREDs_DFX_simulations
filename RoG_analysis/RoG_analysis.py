import mdtraj as md
import pandas as pd
import numpy as np
import glob
import os
import argparse
import sys

def main():
    # Parse command-line arguments
    try:
        parser = argparse.ArgumentParser(description='Compute radius of gyration from MD trajectories and save to CSV.')
        parser.add_argument('--input_folder', required=True, help='Path to folder containing trajectory files')
        parser.add_argument('--output_folder', required=True, help='Path to folder for saving output CSV')
        args = parser.parse_args()
    except argparse.ArgumentError as e:
        print(f"Error parsing command-line arguments: {e}")
        sys.exit(1)

    # Validate input folder
    if not os.path.isdir(args.input_folder):
        print(f"Error: Input folder '{args.input_folder}' does not exist or is not a directory.")
        sys.exit(1)

    # Ensure output folder exists
    try:
        os.makedirs(args.output_folder, exist_ok=True)
    except OSError as e:
        print(f"Error: Could not create output folder '{args.output_folder}': {e}")
        sys.exit(1)

    # Find all production trajectory files
    trajectory_files = sorted(glob.glob(os.path.join(args.input_folder, "production_*.h5")))
    if not trajectory_files:
        print(f"Error: No trajectory files (production_*.h5) found in '{args.input_folder}'.")
        sys.exit(1)

    # Initialize dictionary to store results
    rg_data = {}

    # Process each trajectory file
    for traj_file in trajectory_files:
        print(f"Processing {traj_file}...")
        try:
            # Load trajectory, selecting every 1000th frame
            traj = md.load(traj_file, stride=100)
        except Exception as e:
            print(f"Error: Failed to load trajectory '{traj_file}': {e}")
            continue

        try:
            # Superpose only the selected frames
            traj.superpose(traj, frame=0)

            # Compute radius of gyration
            rg = md.compute_rg(traj)

            # Store results
            replicate_name = os.path.basename(traj_file).replace(".h5", "")
            rg_data[replicate_name] = rg
        except Exception as e:
            print(f"Error: Failed to process trajectory '{traj_file}': {e}")
            continue

    # Check if any data was successfully processed
    if not rg_data:
        print("Error: No trajectories were successfully processed.")
        sys.exit(1)

    # Convert to DataFrame
    try:
        rg_df = pd.DataFrame(rg_data)
    except Exception as e:
        print(f"Error: Failed to create DataFrame: {e}")
        sys.exit(1)

    # Save to CSV
    output_path = os.path.join(args.output_folder, "radius_of_gyration.csv")
    try:
        rg_df.to_csv(output_path, index=False)
        print(f"Analysis complete. Results saved to {output_path}.")
    except Exception as e:
        print(f"Error: Failed to save CSV to '{output_path}': {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
