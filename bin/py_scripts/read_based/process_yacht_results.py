import argparse
import os
import re
import pandas as pd
import shutil

# Set up argument parser
parser = argparse.ArgumentParser(description='Process yacht data.')
parser.add_argument('-d', '--directory', type=str, help='Base folder path', required=True)

# Parse arguments
args = parser.parse_args()

# Assign base folder from command-line argument
# Example directory: /mnt/data1/quantnm/yacht_automation/yacht
base_folder = args.directory

# Folder paths
input_folder = os.path.join(base_folder)
output_folder = os.path.join(base_folder, 'csv')
temp_folder = os.path.join(base_folder, 'temp')

# Create output folders if they don't exist
os.makedirs(output_folder, exist_ok=True)
os.makedirs(temp_folder, exist_ok=True)

# Task 1: Extract sheets and save as CSV
for filename in os.listdir(input_folder):
    if filename.endswith('_yacht.xlsx'):
        xlsx_file = os.path.join(input_folder, filename)
        id = re.search(r'(.*)_yacht\.xlsx', filename).group(1)
        
        excel_data = pd.read_excel(xlsx_file, sheet_name=['min_coverage0.5', 'min_coverage0.1', 'min_coverage0.05'])
        
        for sheet_name, df in excel_data.items():
            csv_filename = f"{id}_{sheet_name}_yacht.csv"
            csv_path = os.path.join(output_folder, csv_filename)
            df.to_csv(csv_path, index=False)
            print(f"Saved {csv_filename}")

# Task 2: Add new columns to the CSV files
for csv_file in os.listdir(output_folder):
    if csv_file.endswith('.csv'):
        csv_path = os.path.join(output_folder, csv_file)
        df = pd.read_csv(csv_path)
        
        # Extract min_coverage<value> using regex
        match = re.search(r'min_coverage\d+\.\d+', csv_file)
        file_name = match.group(0) if match else 'unknown'
        
        # Extract ID (WGS_ID) by removing the min_coverage<value> and the suffix
        wgs_id = csv_file.split('_min_coverage')[0]
        
        df['file_name'] = file_name
        df['WGS_ID'] = wgs_id
        
        new_csv_path = os.path.join(temp_folder, csv_file)
        df.to_csv(new_csv_path, index=False)
        print(f"Saved {new_csv_path}")

csv_files = [pd.read_csv(os.path.join(temp_folder, f)) for f in os.listdir(temp_folder)]
csv_files = [df for df in csv_files if not df.empty and not df.isna().all().all()]
merged_df = pd.concat(csv_files, ignore_index=True)

merged_df.to_csv(os.path.join(base_folder, 'yacht_mergedresults.csv'), index=False, float_format='%.15g')
merged_df.to_excel(os.path.join(base_folder, 'yacht_mergedresults.xlsx'), index=False)

shutil.rmtree(temp_folder)
print(f"Removed {temp_folder}")
