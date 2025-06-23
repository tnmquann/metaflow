import argparse
import os
import pandas as pd

# Set up argument parser
parser = argparse.ArgumentParser(description='Process Sourmash data.')
parser.add_argument('-d', '--directory', type=str, help='Directory path for output', default='.')
parser.add_argument('-i', '--input', type=str, help='Path to the input file need to be processed', required=True)
parser.add_argument('-o', '--output', type=str, help='Name of the output file', required=True)

# Parse arguments
args = parser.parse_args()
input_file = args.input
output_dir = args.directory
output_file_name = args.output

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Function to calculate relative abundance (ref: https://github.com/sourmash-bio/sourmash/issues/461)
def calculate_relative_abundance(row):
    relative_abundance = ((row['unique_intersect_bp']/row['scaled']) * row['median_abund'])/row['total_weighted_hashes']
    return relative_abundance

# Load the merged input file
df = pd.read_csv(input_file)

# Extract WGS_ID from query_name column
df['WGS_ID'] = df['query_name']

# Calculate relative abundance for each unique WGS_ID
df['relative_abundance'] = df.apply(calculate_relative_abundance, axis=1)

# Handle unclassified relative abundance for each unique WGS_ID
unclassified_data = []
for wgs_id in df['WGS_ID'].unique():
    subset = df[df['WGS_ID'] == wgs_id]
    unclassified_relative_abundance = 1 - subset['relative_abundance'].sum()
    unclassified_data.append({
        'name': 'unclassified',
        'WGS_ID': wgs_id,
        'relative_abundance': unclassified_relative_abundance
    })

unclassified_df = pd.DataFrame(unclassified_data)

# Append the unclassified data
df = pd.concat([df, unclassified_df], ignore_index=True)

# Save the processed data
output_file = os.path.join(output_dir, output_file_name)
df.to_csv(output_file, index=False)
print(f"Output file '{output_file}' has been created with the new 'WGS_ID' and 'relative_abundance' columns.")