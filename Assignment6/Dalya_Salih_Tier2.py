import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import NMF
from scipy.cluster.hierarchy import linkage
import sys

def create_symmetric_matrix(file_path):
    # Load the dataset
    data = pd.read_csv(file_path, sep='\t')
    
    # Get all unique protein names to build the symmetric matrix
    proteins = pd.concat([data['protein1'], data['protein2']]).unique()
    protein_to_idx = {protein: idx for idx, protein in enumerate(proteins)}
    n = len(proteins)

    # Initialize the symmetric interaction matrix with zeros
    interaction_matrix = np.zeros((n, n))

    # Fill in the matrix with interaction scores
    for _, row in data.iterrows():
        i, j = protein_to_idx[row['protein1']], protein_to_idx[row['protein2']]
        score = row['interaction_score']
        interaction_matrix[i, j] = score
        interaction_matrix[j, i] = score  # Make it symmetric

    return interaction_matrix, proteins

def apply_nmf_and_cluster(matrix, proteins):
    # Apply NMF with specified parameters
    model = NMF(n_components=10, init='random', random_state=42)
    W = model.fit_transform(matrix)
    H = model.components_

    # Assign each protein to the cluster with the highest value in W
    clusters = W.argmax(axis=1)
    cluster_counts = np.bincount(clusters)
    
    # Find the cluster with the fewest proteins
    target_cluster = np.argmin(cluster_counts)
    
    # Filter proteins in the selected cluster
    selected_indices = np.where(clusters == target_cluster)[0]
    selected_proteins = [proteins[idx] for idx in selected_indices]
    selected_matrix = matrix[np.ix_(selected_indices, selected_indices)]

    return selected_matrix, selected_proteins

def plot_heatmap(selected_matrix, selected_proteins):
    # Perform hierarchical clustering
    linkage_matrix = linkage(selected_matrix, method='ward')
    
    # Create a clustermap with dendrogram
    sns.clustermap(selected_matrix, row_linkage=linkage_matrix, col_linkage=linkage_matrix,
                   xticklabels=selected_proteins, yticklabels=selected_proteins, cmap='viridis')
    
    # Save the heatmap
    plt.savefig("Protein_Cluster_Heatmap.png")
    plt.close()

if __name__ == "__main__":
    # Command-line argument for file path
    file_path = sys.argv[1]

    # Step 1: Load data and create symmetric interaction matrix
    interaction_matrix, proteins = create_symmetric_matrix(file_path)

    # Step 2: Apply NMF and filter proteins by the cluster with the fewest proteins
    selected_matrix, selected_proteins = apply_nmf_and_cluster(interaction_matrix, proteins)

    # Step 3: Plot and save the heatmap with hierarchical clustering
    plot_heatmap(selected_matrix, selected_proteins)
