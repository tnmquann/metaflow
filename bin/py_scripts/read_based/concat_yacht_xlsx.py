import argparse
import os
import pandas as pd

# Set up argument parser
parser = argparse.ArgumentParser(description='Process and merge yacht results.')
parser.add_argument('-d', '--directory', type=str, help='Directory path for input and output', required=True)

# Parse arguments
args = parser.parse_args()

# Use the specified directory path for both input and output
folder_path = args.directory

# Read the yacht_results.csv file
yacht_results = pd.read_csv(os.path.join(folder_path, 'yacht_mergedresults.csv'))

# Extract the required columns
required_cols = ['organism_name', 'p_vals', 'actual_confidence_with_coverage', 'alt_confidence_mut_rate_with_coverage', 'file_name', 'WGS_ID']
yacht_merged = yacht_results[required_cols]

# Convert the required columns to object dtype first, then to string
for col in required_cols:
    yacht_merged.loc[:, col] = yacht_merged.loc[:, col].astype(object).astype(str)

# Merge the records based on organism_name and WGS_ID
yacht_merged = yacht_merged.groupby(['organism_name', 'WGS_ID']).agg(' | '.join).reset_index()

# Save the merged result as yacht_merged_processed.csv
yacht_merged.to_csv(os.path.join(folder_path, 'yacht_merged_processed.csv'), index=False)
print("Saved yacht_merged_processed.csv")