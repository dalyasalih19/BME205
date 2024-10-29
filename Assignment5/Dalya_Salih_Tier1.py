# Firstname_Lastname_Tier1.py
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
import pandas as pd

# Set a random seed for reproducibility
RANDOM_STATE = 42

# Part 1: PCA on MNIST Dataset

# 1. Load the MNIST subset
mnist_data = np.load("MNIST_X_subset.npy")
mnist_labels = np.load("MNIST_y_subset.npy")

# 2. Apply PCA to reduce dimensions to 2
pca_mnist = PCA(n_components=2, random_state=RANDOM_STATE)
mnist_data_2d = pca_mnist.fit_transform(mnist_data)

# 3. Visualize the 2D PCA-transformed data
plt.figure(figsize=(8, 6))
scatter = plt.scatter(mnist_data_2d[:, 0], mnist_data_2d[:, 1], c=mnist_labels, cmap='viridis', alpha=0.5)
plt.colorbar(scatter, label='Digit Label')
plt.title('MNIST PCA (2D)')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.savefig("MNIST_PCA_2D.png")
plt.close()

# 4. Reconstruct an image using the first 2 principal components (index 0)
example_image = mnist_data[0]
example_reduced = pca_mnist.transform([example_image])
example_reconstructed = pca_mnist.inverse_transform(example_reduced)

# Save original and reconstructed images
plt.imshow(example_image.reshape(28, 28), cmap='gray')
plt.title("Original MNIST Image (Index 0)")
plt.savefig("MNIST_original.png")
plt.close()

plt.imshow(example_reconstructed.reshape(28, 28), cmap='gray')
plt.title("Reconstructed MNIST Image (2 PCs)")
plt.savefig("MNIST_reconstructed_2PC.png")
plt.close()

# 5. Reconstruct an image from a chosen 2D point (Digit 1 representation)
chosen_point = np.array([0, -5])  # Example chosen point to resemble "1"
chosen_reconstructed = pca_mnist.inverse_transform([chosen_point])

# Save the reconstructed image of the digit "1"
plt.imshow(chosen_reconstructed.reshape(28, 28), cmap='gray')
plt.title("Reconstructed MNIST '1' from Chosen Coord")
plt.savefig("MNIST_reconstructed_1_from_coord.png")
plt.close()


# Part 2: PCA on Dogs SNP Dataset

# 1. Load the Dogs SNP dataset
dogs_data = np.load("dogs_X.npy")
dogs_clades = np.load("dogs_clades.npy")

# 2. Apply PCA to reduce dimensions to 2
pca_dogs = PCA(n_components=2, random_state=RANDOM_STATE)
dogs_data_2d = pca_dogs.fit_transform(dogs_data)

# 3. Visualize the 2D PCA-transformed data for Dogs SNP
plt.figure(figsize=(8, 6))
scatter = plt.scatter(dogs_data_2d[:, 0], dogs_data_2d[:, 1], c=dogs_clades, cmap='tab10', alpha=0.5)
plt.colorbar(scatter, label='Clade Label')
plt.title('Dogs SNP PCA (2D)')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.savefig("Dogs_PCA_2D.png")
plt.close()


# Part 3: MDS on Molecular Distance Matrix

# 1. Load the molecular distance matrix and atom types
molecule_distances = pd.read_csv("molecule_distances.tsv", sep='\t')
distance_matrix = molecule_distances.iloc[:, 1:].values  # Distance matrix excluding atom label column

# 2. Perform MDS to reconstruct 3D coordinates
mds = MDS(n_components=3, dissimilarity='precomputed', random_state=RANDOM_STATE)
molecule_coords = mds.fit_transform(distance_matrix)

# 3. Output the coordinates to a CSV file
molecule_coords_df = pd.DataFrame(molecule_coords, columns=['X', 'Y', 'Z'])
molecule_coords_df.to_csv("molecule_coordinates.csv", index=False)
