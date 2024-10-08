import os
import csv
import sys
import numpy as np
import statistics

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

def mann_whitney_u_manual(control_data, treatment_data):
    """
    Perform the Mann-Whitney U test manually using ranks, and use NormalDist from statistics to compute p-values.
    :param control_data: Numpy array of normalized control expression (rows = replicates, columns = genes).
    :param treatment_data: Numpy array of normalized treatment expression (rows = replicates, columns = genes).
    :return: List of p-values for each gene.
    """
    num_genes = control_data.shape[1]
    p_values = []

    for gene_idx in range(num_genes):
        # Get the gene expression data for control and treatment
        control_gene_data = control_data[:, gene_idx]
        treatment_gene_data = treatment_data[:, gene_idx]
        
        # Combine the control and treatment data and rank them
        combined_data = np.concatenate([control_gene_data, treatment_gene_data])
        ranks = np.argsort(np.argsort(combined_data)) + 1  # Rank data starting from 1
        
        # Split the ranks back into control and treatment
        ranks_control = ranks[:len(control_gene_data)]
        ranks_treatment = ranks[len(control_gene_data):]
        
        # Compute U-statistic
        U1 = np.sum(ranks_control) - (len(ranks_control) * (len(ranks_control) + 1)) / 2
        U2 = np.sum(ranks_treatment) - (len(ranks_treatment) * (len(ranks_treatment) + 1)) / 2
        U = min(U1, U2)
        
        # Calculate mean and standard deviation for U
        n1 = len(control_gene_data)
        n2 = len(treatment_gene_data)
        mean_U = n1 * n2 / 2
        std_U = np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
        
        # Calculate z-score
        z = (U - mean_U) / std_U
        
        # Calculate two-sided p-value using statistics.NormalDist
        normal_dist = statistics.NormalDist(mu=0, sigma=1)
        p_value = 2 * (1 - normal_dist.cdf(abs(z)))
        p_values.append(p_value)

    return p_values

def main():
    # Get the control and treatment directories from command-line arguments
    if len(sys.argv) != 3:
        print("Usage: python Firstname_Lastname_Tier2.py <control_directory> <treatment_directory>")
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

    # Perform Mann-Whitney U test to calculate p-values (using a manual method)
    p_values = mann_whitney_u_manual(control_data, treatment_data)

    # Collect all data into a list of tuples for sorting
    gene_data = []
    for i, (mc, med_c, mt, med_t, fc, p_val) in enumerate(zip(mean_control, median_control, mean_treatment, median_treatment, log2_fold_change, p_values)):
        gene_data.append((f"Gene{i+1}", mc, med_c, mt, med_t, fc, p_val))

    # Sort gene data by p-value (lowest first)
    gene_data.sort(key=lambda x: x[6])  # Sort by p-value (index 6 in the tuple)

    # Output sorted results
    output_file = 'output_Tier2.txt'
    with open(output_file, 'w') as outfile:
        outfile.write("#gene\tmean_control\tmedian_control\tmean_treatment\tmedian_treatment\tlog2_fold_change\tp_value\n")
        for gene in gene_data:
            outfile.write(f"{gene[0]}\t{gene[1]:.6f}\t{gene[2]:.6f}\t{gene[3]:.6f}\t{gene[4]:.6f}\t{gene[5]:.6f}\t{gene[6]:.6e}\n")

    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    main()
