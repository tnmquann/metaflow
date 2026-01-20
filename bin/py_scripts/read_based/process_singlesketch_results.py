#!/usr/bin/env python3
"""
Process read-based results for a single sample.
This script processes sourmash gather and YACHT results for individual samples.
Unlike batch processing, this doesn't need to concatenate multiple samples.
"""

import argparse
import os
import pandas as pd


def calculate_relative_abundance(row):
    """Calculate relative abundance for sourmash results."""
    relative_abundance = ((row['unique_intersect_bp']/row['scaled']) * row['median_abund'])/row['total_weighted_hashes']
    return relative_abundance


def process_sourmash_results(input_file, output_dir, sample_id):
    """Process sourmash gather results for a single sample."""
    print(f"Processing sourmash results for sample: {sample_id}")
    
    # Load the input file
    df = pd.read_csv(input_file)
    
    if df.empty:
        print(f"Warning: Sourmash results file is empty for sample {sample_id}")
        # Create empty dataframe with required columns
        empty_df = pd.DataFrame(columns=['name', 'WGS_ID', 'relative_abundance'])
        output_file = os.path.join(output_dir, f'{sample_id}_sourmash_processed.csv')
        empty_df.to_csv(output_file, index=False)
        return output_file
    
    # Set WGS_ID as sample_id for single sample processing
    df['WGS_ID'] = sample_id
    
    # Calculate relative abundance
    df['relative_abundance'] = df.apply(calculate_relative_abundance, axis=1)
    
    # Handle unclassified relative abundance
    unclassified_relative_abundance = 1 - df['relative_abundance'].sum()
    unclassified_data = {
        'name': 'unclassified',
        'WGS_ID': sample_id,
        'relative_abundance': unclassified_relative_abundance
    }
    
    unclassified_df = pd.DataFrame([unclassified_data])
    df = pd.concat([df, unclassified_df], ignore_index=True)
    
    # Save the processed data
    output_file = os.path.join(output_dir, f'{sample_id}_sourmash_processed.csv')
    df.to_csv(output_file, index=False)
    print(f"Saved sourmash processed results to: {output_file}")
    
    return output_file


def process_yacht_results(yacht_xlsx_path, output_dir, sample_id):
    """
    Process YACHT results for a single sample.
    Each sample has one YACHT xlsx file with multiple sheets (different min_coverage values).
    No need to concatenate multiple files like in batch processing.
    """
    print(f"Processing YACHT results for sample: {sample_id}")
    
    all_dfs = []
    
    try:
        # Read all sheets from the single YACHT xlsx file
        excel_data = pd.read_excel(yacht_xlsx_path, sheet_name=None)
        
        for sheet_name, df in excel_data.items():
            if df.empty:
                continue
            print(f"  Processing sheet '{sheet_name}'")
            
            # Add metadata columns
            df['file_name'] = sheet_name  # Sheet name contains min_coverage info
            df['WGS_ID'] = sample_id
            
            all_dfs.append(df)
            
    except Exception as e:
        print(f"Error processing YACHT file: {e}")
        return None
    
    if not all_dfs:
        print(f"Warning: No valid YACHT data found for sample {sample_id}")
        return None
    
    # Filter out empty dataframes
    valid_dfs = [df for df in all_dfs if not df.empty and not df.isna().all().all()]
    
    if not valid_dfs:
        print(f"Warning: No valid data found in YACHT sheets for sample {sample_id}")
        return None
    
    # Concatenate all sheets from this single sample's xlsx
    merged_df = pd.concat(valid_dfs, ignore_index=True)
    
    # Process and aggregate YACHT results (combine rows with same organism across sheets)
    required_cols = ['organism_name', 'p_vals', 'actual_confidence_with_coverage', 
                     'alt_confidence_mut_rate_with_coverage', 'file_name', 'WGS_ID']
    
    available_cols = [col for col in required_cols if col in merged_df.columns]
    
    if 'organism_name' not in merged_df.columns or 'WGS_ID' not in merged_df.columns:
        print("Warning: Required columns 'organism_name' or 'WGS_ID' not found in YACHT results")
        return None
    
    yacht_subset = merged_df[available_cols].copy()
    
    # Convert columns to string for aggregation
    for col in available_cols:
        yacht_subset[col] = yacht_subset[col].astype(object).astype(str)
    
    # Aggregate by organism_name and WGS_ID (combine info from different min_coverage sheets)
    yacht_processed = yacht_subset.groupby(['organism_name', 'WGS_ID']).agg(' | '.join).reset_index()
    
    # Save processed results
    yacht_processed_path = os.path.join(output_dir, f'{sample_id}_yacht_processed.csv')
    yacht_processed.to_csv(yacht_processed_path, index=False)
    print(f"Saved YACHT processed results to: {yacht_processed_path}")
    
    return yacht_processed_path


def merge_results(sourmash_csv_path, yacht_csv_path, output_dir, sample_id, 
                  remove_unclassified=False, min_coverage=None):
    """Merge sourmash and YACHT results for a single sample."""
    print(f"Merging results for sample: {sample_id}")
    
    # Load sourmash results
    try:
        sourmash_df = pd.read_csv(sourmash_csv_path)
    except Exception as e:
        print(f"Error loading sourmash results: {e}")
        return
    
    # Rename 'name' column to 'organism_name' for merging
    if 'name' in sourmash_df.columns:
        sourmash_df = sourmash_df.rename(columns={'name': 'organism_name'})
    
    # Load YACHT results if available
    if yacht_csv_path and os.path.exists(yacht_csv_path):
        try:
            yacht_df = pd.read_csv(yacht_csv_path)
            # Merge dataframes
            merged = pd.merge(sourmash_df, yacht_df, on=['organism_name', 'WGS_ID'], how='outer')
        except Exception as e:
            print(f"Warning: Error loading YACHT results: {e}")
            merged = sourmash_df
    else:
        print("Warning: YACHT results not available, using only sourmash results")
        merged = sourmash_df
    
    # Save the initial merged dataframe
    merged_csv_path = os.path.join(output_dir, f'{sample_id}_merged_sourmash_yacht.csv')
    merged_xlsx_path = os.path.join(output_dir, f'{sample_id}_merged_sourmash_yacht.xlsx')
    
    merged.to_csv(merged_csv_path, index=False)
    merged.to_excel(merged_xlsx_path, index=False)
    print(f"Saved merged results to: {merged_csv_path}")
    
    # Define final columns
    final_cols = [
        'organism_name', 'WGS_ID', 'relative_abundance', 'p_vals',
        'actual_confidence_with_coverage', 'alt_confidence_mut_rate_with_coverage', 'file_name'
    ]
    
    # Filter to available columns
    available_final_cols = [col for col in final_cols if col in merged.columns]
    
    # Process with min_coverage filter
    if min_coverage is not None:
        print(f"Processing with --min_coverage {min_coverage}...")
        filtered_df = merged.dropna(subset=['relative_abundance']).copy()
        
        if 'file_name' in filtered_df.columns:
            filtered_df = filtered_df.dropna(subset=['file_name'])
            coverage_str = f"min_coverage{min_coverage}"
            filtered_df = filtered_df[filtered_df['file_name'].astype(str).str.contains(coverage_str, na=False)]
        
        if not filtered_df.empty:
            # Normalize relative abundance
            group_sums = filtered_df.groupby('WGS_ID')['relative_abundance'].transform('sum')
            filtered_df['relative_abundance'] = filtered_df['relative_abundance'] / group_sums
            
            output_df = filtered_df[available_final_cols]
            output_basename = f'{sample_id}_final_mincoverage{min_coverage}'
            
            output_df.to_csv(os.path.join(output_dir, f'{output_basename}.csv'), index=False)
            output_df.to_excel(os.path.join(output_dir, f'{output_basename}.xlsx'), index=False)
            print(f"Saved filtered results: {output_basename}")
        else:
            print("Warning: No data after filtering with min_coverage")
    
    elif remove_unclassified:
        print("Processing with --remove_unclassified...")
        filtered_df = merged.dropna(subset=['relative_abundance']).copy()
        
        if 'file_name' in filtered_df.columns:
            filtered_df = filtered_df.dropna(subset=['file_name'])
        
        if not filtered_df.empty:
            # Normalize relative abundance
            group_sums = filtered_df.groupby('WGS_ID')['relative_abundance'].transform('sum')
            filtered_df['relative_abundance'] = filtered_df['relative_abundance'] / group_sums
            
            output_df = filtered_df[available_final_cols]
            output_basename = f'{sample_id}_final_removeuncl'
            
            output_df.to_csv(os.path.join(output_dir, f'{output_basename}.csv'), index=False)
            output_df.to_excel(os.path.join(output_dir, f'{output_basename}.xlsx'), index=False)
            print(f"Saved filtered results: {output_basename}")
        else:
            print("Warning: No data after filtering")


def main():
    parser = argparse.ArgumentParser(description='Process read-based results for a single sample.')
    parser.add_argument('-s', '--sample_id', type=str, required=True,
                        help='Sample ID')
    parser.add_argument('-g', '--gather_csv', type=str, required=True,
                        help='Path to sourmash gather CSV file')
    parser.add_argument('-y', '--yacht_xlsx', type=str, required=True,
                        help='Path to YACHT xlsx file')
    parser.add_argument('-o', '--output_dir', type=str, required=True,
                        help='Output directory for results')
    parser.add_argument('--remove_unclassified', action='store_true',
                        help='Remove unclassified rows and normalize abundance')
    parser.add_argument('--min_coverage', type=float, default=None,
                        help='Filter by a specific min_coverage value')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Step 1: Process sourmash results
    sourmash_processed_path = process_sourmash_results(
        args.gather_csv, 
        args.output_dir, 
        args.sample_id
    )
    
    # Step 2: Process YACHT results (single file, no concatenation needed)
    yacht_processed_path = process_yacht_results(
        args.yacht_xlsx,
        args.output_dir,
        args.sample_id
    )
    
    # Step 3: Merge results
    merge_results(
        sourmash_processed_path,
        yacht_processed_path,
        args.output_dir,
        args.sample_id,
        remove_unclassified=args.remove_unclassified,
        min_coverage=args.min_coverage
    )
    
    print(f"\nSingle sample processing complete for: {args.sample_id}")


if __name__ == '__main__':
    main()
