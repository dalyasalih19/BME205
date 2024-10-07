import os
import csv
import sys
import numpy as np

def normalize_expression(data):
    """
    Normalize gene expression counts by dividing each gene's count by the total expression in the replicate.
    :param data: List of gene expression data for a single replicate.
    :return: List with normalized gene expression.
    """
    total_expression = sum(data)
    return [x / total_expression for x in data]

def read_and_normalize(directory):
    """
    Read and normalize gene expression data from all CSV files in a given directory.
    Each file contains gene expression values for 10 genes.
    :param directory: Path to the directory containing CSV files.
    :return: Numpy array where each row represents normalized gene expression values for one replicate (10 genes per file).
    """
    all_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.csv')]
    
    if len(all_files) == 0:
        raise FileNotFoundError(f"No CSV files found in directory {directory}.")

    all_gene_data = []

    # Read and normalize each file
    for file in all_files:
        with open(file, mode='r') as infile:
            reader = csv.reader(infile)
            next(reader)  # Skip the header

            gene_expression = []
            for row in reader:  # Process each gene in the file
                gene_expression.append(float(row[1]))  # Extract the gene expression value
            
            # Normalize the expression data for all 10 genes in this file
            normalized = normalize_expression(gene_expression)
            all_gene_data.append(normalized)  # Append the normalized values for this replicate

    # Convert the list of lists to a numpy array
    return np.array(all_gene_data)

def compute_statistics(data):
    """
    Compute the mean and median normalized expression for each gene across all replicates.
    :param data: Numpy array with normalized gene expression for all replicates.
    :return: Tuple of mean and median arrays for each gene.
    """
    mean_expression = np.mean(data, axis=0)
    median_expression = np.median(data, axis=0)
    return mean_expression, median_expression

def calculate_log2_fold_change(mean_control, mean_treatment):
    """
    Calculate the log2 fold change for each gene between treatment and control groups.
    :param mean_control: Numpy array with mean normalized expression for control.
    :param mean_treatment: Numpy array with mean normalized expression for treatment.
    :return: Numpy array with log2 fold change values for each gene.
    """
    return np.log2(mean_treatment / mean_control)

def main():
    # Get the control and treatment directories from command-line arguments
    if len(sys.argv) != 3:
        print("Usage: python Firstname_Lastname_Tier1.py <control_directory> <treatment_directory>")
        sys.exit(1)

    control_directory = sys.argv[1]
    treatment_directory = sys.argv[2]

    # Read and normalize data from both control and treatment groups
    try:
        control_data = read_and_normalize(control_directory)
        treatment_data = read_and_normalize(treatment_directory)
    except (FileNotFoundError, ValueError) as e:
        print(f"Error reading files: {e}")
        sys.exit(1)

    # Compute mean and median normalized expression for each gene
    mean_control, median_control = compute_statistics(control_data)
    mean_treatment, median_treatment = compute_statistics(treatment_data)

    # Compute log2 fold change for each gene
    log2_fold_change = calculate_log2_fold_change(mean_control, mean_treatment)

    # Output results
    output_file = 'output_Tier1.txt'
    with open(output_file, 'w') as outfile:
        outfile.write("#gene\tmean_control\tmedian_control\tmean_treatment\tmedian_treatment\tlog2_fold_change\n")
        for i, (mc, med_c, mt, med_t, fc) in enumerate(zip(mean_control, median_control, mean_treatment, median_treatment, log2_fold_change)):
            outfile.write(f"Gene{i+1}\t{mc:.6f}\t{med_c:.6f}\t{mt:.6f}\t{med_t:.6f}\t{fc:.6f}\n")

    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    main()
