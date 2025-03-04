import argparse
import os
import pandas as pd

# Set up argument parser
parser = argparse.ArgumentParser(description='Merge Sourmash and Yacht results.')
parser.add_argument('-d1', '--directory', type=str, help='Output directory for the merged files', required=True)
parser.add_argument('-d2', '--sourmash_dir', type=str, help='Directory where Sourmash CSV files are stored', required=True)
parser.add_argument('-d3', '--yacht_dir', type=str, help='Directory where Yacht CSV files are stored', required=True)

# Parse arguments
args = parser.parse_args()

# Use the specified directory paths
directory = args.directory
sourmash_directory = args.sourmash_dir
yacht_directory = args.yacht_dir

# Load the CSV files
sourmash_csv = pd.read_csv(os.path.join(sourmash_directory, 'sourmash_mergedresults.csv'))
yacht_csv = pd.read_csv(os.path.join(yacht_directory, 'yacht_merged_processed.csv'))

# Rename the 'name' column in 'sourmash_benchmarking_processed.csv' to 'organism_name'
sourmash_csv = sourmash_csv.rename(columns={'name': 'organism_name'})

# Merge the dataframes based on 'organism_name' and 'WGS_ID'
merged = pd.merge(sourmash_csv, yacht_csv, on=['organism_name', 'WGS_ID'], how='outer')

# Save the merged dataframe to a new CSV file and an Excel file
merged.to_csv(os.path.join(directory, 'merged_sourmash_yacht.csv'), index=False)
merged.to_excel(os.path.join(directory, 'merged_sourmash_yacht.xlsx'), index=False)