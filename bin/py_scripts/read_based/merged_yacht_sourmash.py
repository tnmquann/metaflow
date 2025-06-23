import argparse
import os
import pandas as pd

# --- Functions ---

def process_and_save(df, output_basename, output_dir, final_cols):
    """
    Normalizes the relative_abundance column, selects final columns, and saves the result.

    Args:
        df (pd.DataFrame): The input, filtered DataFrame.
        output_basename (str): The base name for the output files (without extension).
        output_dir (str): The directory to save the output files.
        final_cols (list): The list of columns to keep in the final output.
    """
    # Check if the DataFrame is empty after filtering
    if df.empty:
        print(f"Warning: DataFrame is empty after filtering. No output file '{output_basename}' will be created.")
        return

    # (2) Normalize the relative_abundance column for each WGS_ID
    # Calculate the sum of 'relative_abundance' for each 'WGS_ID' group
    group_sums = df.groupby('WGS_ID')['relative_abundance'].transform('sum')
    # Normalize by dividing by the group sum
    df['relative_abundance'] = df['relative_abundance'] / group_sums

    # (3) Keep only the required columns
    final_df = df[final_cols]

    # (4) Export the result files
    csv_path = os.path.join(output_dir, f"{output_basename}.csv")
    xlsx_path = os.path.join(output_dir, f"{output_basename}.xlsx")

    final_df.to_csv(csv_path, index=False)
    final_df.to_excel(xlsx_path, index=False)
    print(f"Generated filtered file: {csv_path}")
    print(f"Generated filtered file: {xlsx_path}")


# --- Main Script ---

# Set up argument parser
parser = argparse.ArgumentParser(description='Merge Sourmash and Yacht results.')
parser.add_argument('-d1', '--directory', type=str, help='Output directory for the merged files', required=True)
parser.add_argument('-d2', '--sourmash_dir', type=str, help='Directory where Sourmash CSV files are stored', required=True)
parser.add_argument('-d3', '--yacht_dir', type=str, help='Directory where Yacht CSV files are stored', required=True)
parser.add_argument('--remove_unclassified', action='store_true', help='Remove unclassified rows and normalize abundance.')
parser.add_argument('--min_coverage', type=float, help='Filter by a specific min_coverage value, remove unclassified, and normalize.')

# Parse arguments
args = parser.parse_args()

# Use the specified directory paths
directory = args.directory
sourmash_directory = args.sourmash_dir
yacht_directory = args.yacht_dir

# Load the CSV files
try:
    sourmash_csv = pd.read_csv(os.path.join(sourmash_directory, 'sourmash_mergedresults.csv'))
    yacht_csv = pd.read_csv(os.path.join(yacht_directory, 'yacht_merged_processed.csv'))
except FileNotFoundError as e:
    print(f"Error: Could not find input file. {e}")
    exit(1)

# Rename the 'name' column in sourmash_csv to 'organism_name' for merging
if 'name' in sourmash_csv.columns:
    sourmash_csv = sourmash_csv.rename(columns={'name': 'organism_name'})

# Merge the dataframes based on 'organism_name' and 'WGS_ID'
merged = pd.merge(sourmash_csv, yacht_csv, on=['organism_name', 'WGS_ID'], how='outer')

# Save the initial merged dataframe
merged.to_csv(os.path.join(directory, 'merged_sourmash_yacht.csv'), index=False)
merged.to_excel(os.path.join(directory, 'merged_sourmash_yacht.xlsx'), index=False)
print(f"Saved initial merged file to {os.path.join(directory, 'merged_sourmash_yacht.csv')}")


# Define the final columns to keep for the processed files
FINAL_COLS = [
    'organism_name', 'WGS_ID', 'relative_abundance', 'p_vals',
    'actual_confidence_with_coverage', 'alt_confidence_mut_rate_with_coverage', 'file_name'
]

# Process the logic for the additional options
# The --min_coverage option is more specific and includes the filtering of --remove_unclassified
if args.min_coverage is not None:
    print(f"Processing with --min_coverage {args.min_coverage}...")
    # (1) Remove rows where 'relative_abundance' or 'file_name' are empty
    filtered_df = merged.dropna(subset=['relative_abundance', 'file_name']).copy()

    # (2) Keep only rows containing the specified min_coverage string
    coverage_str = f"min_coverage{args.min_coverage}"
    # Ensure 'file_name' column is string type to use .str methods
    filtered_df = filtered_df[filtered_df['file_name'].astype(str).str.contains(coverage_str, na=False)]

    # (3 & 4) Normalize and save
    output_basename = f"merged_sourmash_yacht_removeuncl_mincoverage{args.min_coverage}"
    process_and_save(filtered_df, output_basename, directory, FINAL_COLS)

elif args.remove_unclassified:
    print("Processing with --remove_unclassified...")
    # (1) Remove rows where 'relative_abundance' or 'file_name' are empty
    filtered_df = merged.dropna(subset=['relative_abundance', 'file_name']).copy()

    # (2 & 3 & 4) Normalize, select columns, and save
    output_basename = "merged_sourmash_yacht_removeuncl"
    process_and_save(filtered_df, output_basename, directory, FINAL_COLS)