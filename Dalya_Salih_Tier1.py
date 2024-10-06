import os
import pandas as pd
import numpy as np

def normalize_expression(df):
    """
    Normalize gene expression counts by dividing each gene's count by the total expression in the replicate.
    :param df: DataFrame containing gene expression data for a single replicate.
    :return: DataFrame with normalized gene expression.
    """
    return df.div(df.sum(axis=1), axis=0)

def read_and_normalize(directory):
    """
    Read and normalize gene expression data from all files in a given directory.
    :param directory: Path to the directory containing CSV files.
    :return: DataFrame where each row represents a replicate's normalized gene expression.
    """
    all_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.csv')]
    expression_data = []
    
    for file in all_files:
        df = pd.read_csv(file, index_col=0)  # Assumes gene IDs are in the first column
        normalized_df = normalize_expression(df)
        expression_data.append(normalized_df)
    
    return pd.concat(expression_data)

def compute_statistics(normalized_data):
    """
    Compute the mean and median normalized expression for each gene.
    :param normalized_data: DataFrame with normalized gene expression for multiple replicates.
    :return: DataFrame with mean and median for each gene.
    """
    mean_expression = normalized_data.mean(axis=0)
    median_expression = normalized_data.median(axis=0)
    return mean_expression, median_expression

def calculate_log2_fold_change(mean_control, mean_treatment):
    """
    Calculate the log2 fold change for each gene between treatment and control groups.
    :param mean_control: Series with mean normalized expression for control.
    :param mean_treatment: Series with mean normalized expression for treatment.
    :return: Series with log2 fold change values for each gene.
    """
    return np.log2(mean_treatment / mean_control)

def main():
    # Step 1: Read and normalize data from both control and treatment groups
    control_directory = 'control_files'
    treatment_directory = 'treatment_files'
    
    control_data = read_and_normalize(control_directory)
    treatment_data = read_and_normalize(treatment_directory)
    
    # Step 2: Compute mean and median normalized expression for each gene
    mean_control, median_control = compute_statistics(control_data)
    mean_treatment, median_treatment = compute_statistics(treatment_data)
    
    # Step 3: Compute log2 fold change for each gene
    log2_fold_change = calculate_log2_fold_change(mean_control, mean_treatment)
    
    # Step 4: Combine the results into a single DataFrame
    result = pd.DataFrame({
        'gene': mean_control.index,
        'mean_control': mean_control.values,
        'median_control': median_control.values,
        'mean_treatment': mean_treatment.values,
        'median_treatment': median_treatment.values,
        'log2_fold_change': log2_fold_change.values
    })
    
    # Step 5: Sort by log2 fold change (lowest to highest)
    result = result.sort_values(by='log2_fold_change')
    
    # Step 6: Write the results to a tab-delimited file
    output_file = 'gene_expression_log2_fold_change.tsv'
    result.to_csv(output_file, sep='\t', index=False, header=[
        'gene', 'mean_control', 'median_control', 'mean_treatment', 'median_treatment', 'log2_fold_change'
    ])
    
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    main()
