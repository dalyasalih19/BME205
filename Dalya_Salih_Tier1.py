import os
import sys
import pandas as pd
import numpy as np

def normalize_expression(df):
    """
    Normalize gene expression counts for a single replicate (DataFrame).
    This is done by dividing each gene's expression count by the total expression of all genes in that replicate.
    
    :param df: DataFrame containing gene expression data for one replicate.
    :return: DataFrame with normalized gene expression values.
    """
    return df.div(df.sum(axis=1), axis=0)

def read_and_normalize(directory):
    """
    Read and normalize gene expression data from all CSV files in a given directory.
    Each file represents an experimental unit with expression data for 10 genes.
    
    :param directory: Path to the directory containing CSV files.
    :return: DataFrame where each row represents normalized gene expression for a replicate.
    """
    # Check if the directory exists
    if not os.path.exists(directory):
        raise FileNotFoundError(f"Directory {directory} does not exist.")
    
    # Get all CSV files in the directory
    all_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.csv')]
    
    if len(all_files) == 0:
        raise FileNotFoundError(f"No CSV files found in directory {directory}.")

    expression_data = []
    
    # Read and normalize each file
    for file in all_files:
        try:
            df = pd.read_csv(file, index_col=0)  # Assumes gene IDs are in the first column
            normalized_df = normalize_expression(df)
            expression_data.append(normalized_df)
        except Exception as e:
            raise ValueError(f"Error processing file {file}: {e}")
    
    # Concatenate data from all replicates into one DataFrame
    return pd.concat(expression_data)

def compute_statistics(normalized_data):
    """
    Compute the mean and median normalized expression for each gene across all replicates.
    
    :param normalized_data: DataFrame containing normalized gene expression data for all replicates.
    :return: Two Series: mean and median normalized expression values for each gene.
    """
    mean_expression = normalized_data.mean(axis=0)
    median_expression = normalized_data.median(axis=0)
    return mean_expression, median_expression

def calculate_log2_fold_change(mean_control, mean_treatment):
    """
    Calculate the log2 fold change for each gene between treatment and control groups.
    
    :param mean_control: Series containing the mean normalized expression for control.
    :param mean_treatment: Series containing the mean normalized expression for treatment.
    :return: Series with log2 fold change values for each gene.
    """
    # Calculate fold change directly without adding a small constant
    fold_change = mean_treatment / mean_control
    
    # Return log2 of the fold change
    return np.log2(fold_change)

def main():
    """
    Main function to process gene expression data from control and treatment directories.
    It normalizes the data, computes statistics, calculates log2 fold change, and outputs the results.
    """
    # Ensure that the script receives two arguments for control and treatment directories
    if len(sys.argv) != 3:
        print("Usage: python your_assignment.py <control_directory> <treatment_directory>")
        sys.exit(1)

    control_directory = sys.argv[1]
    treatment_directory = sys.argv[2]

    # Step 1: Read and normalize data from both control and treatment groups
    try:
        control_data = read_and_normalize(control_directory)
        treatment_data = read_and_normalize(treatment_directory)
    except (FileNotFoundError, ValueError) as e:
        print(f"Error reading files: {e}")
        sys.exit(1)

    # Step 2: Compute mean and median normalized expression for each gene
    mean_control, median_control = compute_statistics(control_data)
    mean_treatment, median_treatment = compute_statistics(treatment_data)

    # Step 3: Compute log2 fold change for each gene
    log2_fold_change = calculate_log2_fold_change(mean_control, mean_treatment)

    # Step 4: Combine the results into a single DataFrame for output
    result = pd.DataFrame({
        'gene': mean_control.index,  # Gene IDs
        'mean_control': mean_control.values,  # Mean control expression
        'median_control': median_control.values,  # Median control expression
        'mean_treatment': mean_treatment.values,  # Mean treatment expression
        'median_treatment': median_treatment.values,  # Median treatment expression
        'log2_fold_change': log2_fold_change.values  # Log2 fold change
    })

    # Step 5: Sort by log2 fold change (lowest to highest)
    result = result.sort_values(by='log2_fold_change')

    # Step 6: Write the results to a tab-delimited file
    output_file = 'gene_expression_log2_fold_change.tsv'
    try:
        result.to_csv(output_file, sep='\t', index=False, header=[
            'gene', 'mean_control', 'median_control', 'mean_treatment', 'median_treatment', 'log2_fold_change'
        ])
        print(f"Results saved to {output_file}")
    except Exception as e:
        print(f"Error writing to file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

