import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import NMF
import sys

def apply_nmf(dog_X_path, dog_clades_path):
    # Load dataset
    dogs_X = np.load(dog_X_path)
    dog_clades = np.load(dog_clades_path, allow_pickle=True)

    # Apply NMF with specified parameters
    model = NMF(n_components=5, init='random', random_state=42)
    W = model.fit_transform(dogs_X)
    H = model.components_

    # Normalize W to represent proportions
    W_normalized = W / W.sum(axis=1, keepdims=True)

    # Determine the dominant cluster for each dog
    dominant_clusters = W_normalized.argmax(axis=1)
    dominant_proportions = W_normalized.max(axis=1)

    # Sort indices based on dominant cluster and proportion
    sorted_indices = np.lexsort((dominant_proportions, dominant_clusters))
    W_sorted = W_normalized[sorted_indices]

    # Create a stacked bar plot
    plt.figure(figsize=(14, 7))
    plt.stackplot(range(W_sorted.shape[0]), W_sorted.T, labels=[f"Cluster {i+1}" for i in range(W_sorted.shape[1])])
    plt.xlabel("Dogs")
    plt.ylabel("Proportion")
    plt.title("NMF Components Proportion for Dogs")
    plt.legend(loc='upper right')
    plt.savefig("NMF_Dogs.png")
    plt.close()

    # Identify dominant component for Basenji samples
    basenji_indices = np.where(dog_clades == "Basenji")[0]
    basenji_dominant_clusters = dominant_clusters[basenji_indices]
    basenji_dominant = np.bincount(basenji_dominant_clusters).argmax()
    print(f"Dominant component in Basenji samples: Cluster {basenji_dominant + 1}")

if __name__ == "__main__":
    dog_X_path = sys.argv[1]
    dog_clades_path = sys.argv[2]
    apply_nmf(dog_X_path, dog_clades_path)

