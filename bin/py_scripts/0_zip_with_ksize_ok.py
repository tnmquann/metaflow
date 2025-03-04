import pandas as pd
import os
import zipfile
import shutil
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description='Process manysketch output')
parser.add_argument('-d', '--directory', type=str, help='Directory path for extracted manysketch output', required=True)
parser.add_argument('-m', '--manifest_dir', type=str, help='Name of the manifest file (without .csv extension)', default='SOURMASH-MANIFEST')

# Parse arguments
args = parser.parse_args()

# Use the specified directory path for both input and output
base_dir = args.directory
manifest_file_name = args.manifest_dir + '.csv'

# Load the manifest file, skipping the first line
manifest_path = os.path.join(base_dir, manifest_file_name)
manifest_gather = pd.read_csv(manifest_path, skiprows=1)

# Create a new column "name_with_ksize"
manifest_gather['name_with_ksize'] = manifest_gather['name'] + '_' + manifest_gather['ksize'].astype(str)

# Function to create zip files with the specified contents
def create_zip_file(row, base_dir):
    internal_location = row['internal_location']
    name_with_ksize = row['name_with_ksize']
    
    # Paths
    file_path = os.path.join(base_dir, internal_location)
    sig_zip_folder = os.path.join(base_dir, 'zip_files')
    zip_file_path = os.path.join(sig_zip_folder, f"{name_with_ksize}.sig.zip")
    
    # Ensure the sig_zip folder exists
    os.makedirs(sig_zip_folder, exist_ok=True)
    
    # Remove "signatures/" from internal_location
    internal_location_short = internal_location.replace("signatures/", "")
    
    # Create manifest.csv content
    manifest_content = (
        "# SOURMASH-MANIFEST-VERSION: 1.0\n"
        "internal_location,md5,md5short,ksize,moltype,num,scaled,n_hashes,with_abundance,name,filename\n"
        f"{row['internal_location']},{row['md5']},{row['md5short']},{row['ksize']},{row['moltype']},{row['num']},{row['scaled']},{row['n_hashes']},{row['with_abundance']},{row['name']},{name_with_ksize}\n"
    )
    
    # Create a temporary directory to hold the files to be zipped
    temp_dir = os.path.join(sig_zip_folder, 'temp')
    os.makedirs(temp_dir, exist_ok=True)
    
    # Create a 'signatures' folder inside the temporary directory
    sigg_folder = os.path.join(temp_dir, 'signatures')
    os.makedirs(sigg_folder, exist_ok=True)
    
    # Move the internal_location file to the 'signatures' folder
    shutil.copy(file_path, sigg_folder)
    
    # Write manifest.csv to the temp directory
    with open(os.path.join(temp_dir, manifest_file_name), 'w') as f:
        f.write(manifest_content)
    
    # Create the zip file
    with zipfile.ZipFile(zip_file_path, 'w') as zipf:
        # Add the internal_location file from the 'signatures' folder
        zipf.write(os.path.join(sigg_folder, os.path.basename(internal_location)), f"signatures/{os.path.basename(internal_location)}")
        # Add the manifest.csv file
        zipf.write(os.path.join(temp_dir, manifest_file_name), manifest_file_name)
    
    # Clean up the temporary directory
    shutil.rmtree(temp_dir)

# Loop over each row in manifest_gather and create the zip files
manifest_gather.apply(create_zip_file, axis=1, base_dir=base_dir)
