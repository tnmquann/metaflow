import argparse
import os
import re
import pandas as pd
import shutil

# Set up argument parser to accept command-line arguments
parser = argparse.ArgumentParser(description='Process YACHT data from Excel files.')
parser.add_argument('-d', '--directory', type=str, help='Base directory path for input and output', required=True)

# Parse the command-line arguments provided by the user
args = parser.parse_args()

# Assign base folder from the command-line argument
base_folder = args.directory
input_folder = os.path.join(base_folder)

# Create an output directory for individual CSV files from each sheet
csv_output_folder = os.path.join(base_folder, 'csv')
os.makedirs(csv_output_folder, exist_ok=True)

# --- Main Logic ---
# This list will hold all the dataframes from all sheets of all Excel files
all_dfs = []

# Iterate over each file in the specified input directory
for filename in os.listdir(input_folder):
    # Process only files that end with '_yacht.xlsx'
    if filename.endswith('_yacht.xlsx'):
        xlsx_file = os.path.join(input_folder, filename)
        # Extract the WGS_ID from the filename using a regular expression
        wgs_id_match = re.search(r'(.*)_yacht\.xlsx', filename)
        if not wgs_id_match:
            print(f"Warning: Could not extract WGS_ID from filename {filename}. Skipping.")
            continue
        wgs_id = wgs_id_match.group(1)
        
        try:
            # Read all sheets from the Excel file into a dictionary of dataframes
            # Setting sheet_name=None tells pandas to read all sheets
            excel_data = pd.read_excel(xlsx_file, sheet_name=None)
            
            # Iterate over each sheet (key-value pair) in the loaded Excel data
            for sheet_name, df in excel_data.items():
                print(f"Processing sheet '{sheet_name}' from {filename}")
                
                # Add metadata columns to the dataframe
                df['file_name'] = sheet_name  # The name of the source sheet
                df['WGS_ID'] = wgs_id         # The WGS_ID extracted from the filename
                
                # Save the current sheet as a separate CSV file
                # Create a safe filename by replacing characters that are not letters, numbers, or ._-
                safe_sheet_name = re.sub(r'[^a-zA-Z0-9_.-]', '_', sheet_name)
                individual_csv_path = os.path.join(csv_output_folder, f"{wgs_id}_{safe_sheet_name}.csv")
                df.to_csv(individual_csv_path, index=False, float_format='%.15g')
                print(f"Saved individual sheet to {individual_csv_path}")

                # Append the processed dataframe to the list for final merging
                all_dfs.append(df)

        except Exception as e:
            # Print an error message if a file cannot be processed
            print(f"Error: Could not process file {filename}. Reason: {e}")

# --- Final Merging and Saving ---
# Proceed only if there are dataframes to process
if all_dfs:
    # Filter out any dataframes that are completely empty or contain only NA values
    valid_dfs = [df for df in all_dfs if not df.empty and not df.isna().all().all()]
    
    # Proceed only if there are valid dataframes left after filtering
    if valid_dfs:
        # Concatenate all valid dataframes into a single large dataframe
        merged_df = pd.concat(valid_dfs, ignore_index=True)

        # Define the output paths for the final merged files
        output_csv_path = os.path.join(base_folder, 'yacht_mergedresults.csv')
        output_xlsx_path = os.path.join(base_folder, 'yacht_mergedresults.xlsx')
        
        # Save the merged dataframe to both CSV and Excel formats
        merged_df.to_csv(output_csv_path, index=False, float_format='%.15g')
        merged_df.to_excel(output_xlsx_path, index=False)
        
        print(f"Successfully saved merged results to {output_csv_path} and {output_xlsx_path}")
    else:
        # Message if no valid data was found in any of the sheets
        print("Warning: No valid data found in any sheets to merge.")
else:
    # Message if no Excel files were found or no dataframes were created
    print("Warning: No '_yacht.xlsx' files found or no dataframes to process.")

# --- Cleanup ---
# Clean up other temporary directories if they exist (e.g., from previous script versions)
# The 'csv' folder is intentionally kept as part of the desired output.
temp_folder = os.path.join(base_folder, 'temp')
if os.path.exists(temp_folder):
    shutil.rmtree(temp_folder)
    print(f"Removed temporary directory: {temp_folder}")
