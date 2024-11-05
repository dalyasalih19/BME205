import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import NMF
import sys

def apply_nmf(dog_X_path, dog_clades_path):
    # Load the Dogs dataset
    dogs_X = np.load(dog_X_path)
    dog_clades = np.load(dog_clades_path, allow_pickle=True)

    # Apply NMF with specified parameters
    model = NMF(n_components=5, init='random', random_state=42)
    W = model.fit_transform(dogs_X)
    H = model.components_

    # Normalize W to make each row sum to 1 (representing proportions)
    W_normalized = W / W.sum(axis=1, keepdims=True)

    # Determine the dominant (largest value) cluster for each dog
    dominant_clusters = W_normalized.argmax(axis=1)
    dominant_proportions = W_normalized.max(axis=1)

    # Sort indices based on dominant cluster and proportion for plotting
    sorted_indices = np.lexsort((dominant_proportions, dominant_clusters))
    W_sorted = W_normalized[sorted_indices]

    # Create a stacked plot for the proportions of each component
    plt.figure(figsize=(14, 7))
    plt.stackplot(range(W_sorted.shape[0]), W_sorted.T, labels=[f"Cluster {i+1}" for i in range(W_sorted.shape[1])])
    plt.xlabel("Dogs")
    plt.ylabel("Proportion")
    plt.title("NMF Components Proportion for Dogs")
    plt.legend(loc='upper right')
    plt.savefig("NMF_Dogs.png")
    plt.close()

    # Identify the dominant component for all Basenji samples
    basenji_indices = np.where(dog_clades == "Basenji")[0]
    if basenji_indices.size > 0:
        basenji_dominant_clusters = dominant_clusters[basenji_indices]
        basenji_dominant = np.bincount(basenji_dominant_clusters).argmax()
        print(basenji_dominant + 1)  # Print the dominant component for Basenji samples
    else:
        print("No Basenji samples found.")

if __name__ == "__main__":
    # Ensure we get file paths from command-line arguments
    dog_X_path = sys.argv[1]
    dog_clades_path = sys.argv[2]
    apply_nmf(dog_X_path, dog_clades_path)
