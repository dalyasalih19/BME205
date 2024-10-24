import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from collections import Counter

# Load the MNIST data from the provided .npy files -- images stored in a 1D, flattened format
MNIST_X = np.load('MNIST_X_subset.npy')
MNIST_y = np.load('MNIST_y_subset.npy')

# Reshape and visualize the first image (index 0) -- images are originally 28x28 pixels, so we reshape them back to 2D to display
first_image = MNIST_X[0].reshape(28, 28)
plt.imshow(first_image, cmap='gray')  # Display the image in grayscale
plt.title(f"Label: {MNIST_y[0]}")  # Add the corresponding label as the title
plt.show()  # Show the plot

# Function to visualize and save the centroids as images
# Always use 1 row, and columns equal to K
def visualize_centroids(centroids, k, filename):
    rows = 1  # Always set the number of rows to 1
    cols = k  # Set the number of columns equal to K (the number of clusters)
    
    for i, centroid in enumerate(centroids):
        # Create a subplot for each centroid (i+1 to avoid indexing from 0)
        plt.subplot(rows, cols, i + 1)
        # Reshape each centroid back into a 28x28 image and display it
        plt.imshow(centroid.reshape(28, 28), cmap='gray')
        plt.axis('off')  # Turn off the axes for cleaner visualization
    plt.suptitle(f'Centroids for K={k}')  # Add a title with the current value of K
    plt.savefig(filename)  # Save the plot to a file
    plt.show()  # Display the plot

# Perform K-Means clustering with K=10 -- KMeans algorithm tries to partition the dataset into 10 clusters (one for each digit)
kmeans_k10 = KMeans(n_clusters=10, random_state=0)
kmeans_k10.fit(MNIST_X)  # Fit the model to the dataset (this performs the clustering)
centroids_k10 = kmeans_k10.cluster_centers_  # Get the centroids of the clusters
# Visualize the centroids for K=10 in a 1x10 grid and save the image as 'centroids_k10.png'
visualize_centroids(centroids_k10, 10, 'centroids_k10.png')

# Perform K-Means clustering with K=11
# Here we try to partition the dataset into 11 clusters, an additional cluster compared to K=10
kmeans_k11 = KMeans(n_clusters=11, random_state=0)
kmeans_k11.fit(MNIST_X)  # Fit the model to the dataset
centroids_k11 = kmeans_k11.cluster_centers_  # Get the centroids of the clusters
# Visualize the centroids for K=11 in a 1x11 grid and save the image as 'centroids_k11.png'
visualize_centroids(centroids_k11, 11, 'centroids_k11.png')

# Function to calculate clustering error
# - kmeans: the fitted KMeans model
# - labels: the true labels for the dataset
def clustering_error(kmeans, labels):
    error = 0  # Initialize the total error counter
    for i in range(kmeans.n_clusters):  # Loop over each cluster
        # Get the indices of the samples assigned to the current cluster (i)
        cluster_indices = np.where(kmeans.labels_ == i)[0]
        # Get the true labels of those samples
        cluster_labels = labels[cluster_indices]
        # Find the majority label in the current cluster (most common digit)
        majority_label = Counter(cluster_labels).most_common(1)[0][0]
        # Count how many samples in the cluster are not equal to the majority label
        error += np.sum(cluster_labels != majority_label)
    return error  # Return the total number of misclassified samples

# Compute and print the clustering error for K=10
error_k10 = clustering_error(kmeans_k10, MNIST_y)
# Print the error for K=10 in the required format
print(f"K=10 Error={error_k10}")

# Compute and print the clustering error for K=11
error_k11 = clustering_error(kmeans_k11, MNIST_y)
# Print the error for K=11 in the required format
print(f"K=11 Error={error_k11}")
