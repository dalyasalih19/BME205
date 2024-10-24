import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import pairwise_distances
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from collections import Counter

# Load the MNIST data
MNIST_X = np.load('MNIST_X_subset.npy')
MNIST_y = np.load('MNIST_y_subset.npy')

# Perform hierarchical clustering with average linkage
# We use 'average' linkage and compute the linkage matrix
Z = linkage(MNIST_X, method='average', metric='euclidean')

# Create the dendrogram and save the plot
plt.figure(figsize=(10, 7))
dn = dendrogram(Z, truncate_mode='lastp', p=10, show_leaf_counts=True)

# Label the clusters with the most common ground truth class digit
cluster_labels = fcluster(Z, t=10, criterion='maxclust')  # Create flat clusters for k=10
# Add cluster labels for the 10 terminal nodes
for i in range(1, 11):
    # Get the indices of the samples in the current cluster
    cluster_indices = np.where(cluster_labels == i)[0]
    # Get the actual labels of these samples
    cluster_digits = MNIST_y[cluster_indices]
    # Find the most common label (majority class) in this cluster
    majority_digit = Counter(cluster_digits).most_common(1)[0][0]
    # Label the cluster in the dendrogram
    plt.text(i*10, 0, f'Majority: {majority_digit}', rotation=90, fontsize=10)

plt.title('MNIST Hierarchical Clustering Dendrogram (k=10)')
plt.savefig('MNIST_dendrogram.png')
plt.show()

# Function to calculate clustering error (similar to Tier 1)
def clustering_error(cluster_labels, true_labels):
    error = 0
    for i in range(1, 11):  # Loop over each cluster
        cluster_indices = np.where(cluster_labels == i)[0]
        cluster_true_labels = true_labels[cluster_indices]
        majority_label = Counter(cluster_true_labels).most_common(1)[0][0]
        error += np.sum(cluster_true_labels != majority_label)
    return error

# Calculate clustering error for hierarchical clustering with k=10
error_hierarchical = clustering_error(cluster_labels, MNIST_y)
print(f"Hierarchical Clustering Error: {error_hierarchical}")

# Save the explanation to a file
with open('MNIST_paragraph.txt', 'w') as f:
    explanation = (
        "The hierarchical clustering groups similar digits together, but the clusters "
        "don't perfectly align with the actual digit labels. For example, some clusters "
        "may contain a mix of digits that have similar shapes (like 4 and 9). This partially "
        "matches expectations, as digits with similar features may end up in the same cluster."
    )
    f.write(explanation)
