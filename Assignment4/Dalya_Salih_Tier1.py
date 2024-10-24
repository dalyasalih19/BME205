import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from collections import Counter

# We load the subset of the MNIST dataset, which contains 6,000 images represented as 1D arrays and their corresponding digit labels
MNIST_X = np.load('MNIST_X_subset.npy')
MNIST_y = np.load('MNIST_y_subset.npy')

# We reshape the first image back to 28x28 pixels to visualize it and turn off the axes for a cleaner view.
first_image = MNIST_X[0].reshape(28, 28)
plt.imshow(first_image, cmap='gray')
plt.axis('off')  # Remove axes for the first image
plt.title(f"Label: {MNIST_y[0]}")  # Add the label as the title
plt.show()  # Show the plot

# This function visualizes the centroids of the clusters generated by K-Means. Each centroid is reshaped back into a 28x28 image, and we save the visualization to an image file.
def visualize_centroids(centroids, k, filename):
    rows = 1  # Always 1 row
    cols = k  # Number of columns equal to number of clusters
    
    # Set figure size to 12 inches wide and 4 inches tall
    plt.figure(figsize=(12, 4)) 
    
    for i, centroid in enumerate(centroids):
        # Create a subplot for each centroid
        ax = plt.subplot(rows, cols, i + 1)
        # Reshape each centroid back into a 28x28 image and display it
        ax.imshow(centroid.reshape(28, 28), cmap='gray')
        ax.axis('off')  # Turn off the axes for cleaner visualization
        ax.set_title(f'Centroid {i + 1}')  # Add title for each individual centroid

    plt.suptitle(f'Centroids for K={k}')  # Add a title with the current value of K at the top
    plt.savefig(filename)  # Save the plot to a file
    plt.show()  # Display the plot

# We perform K-Means clustering with 10 clusters, each representing one digit. We then visualize the centroids of these clusters as 28x28 images.
kmeans_k10 = KMeans(n_clusters=10, random_state=0)
kmeans_k10.fit(MNIST_X)  # Fit the model to the dataset (this performs the clustering)
centroids_k10 = kmeans_k10.cluster_centers_  # Get the centroids of the clusters
# Visualize the centroids for K=10 in a 1x10 grid and save the image as 'centroids_k10.png'
visualize_centroids(centroids_k10, 10, 'centroids_k10.png')

# We run K-Means with 11 clusters to explore how adding one extra cluster affects the centroids and their visual representation.
kmeans_k11 = KMeans(n_clusters=11, random_state=0)
kmeans_k11.fit(MNIST_X)  # Fit the model to the dataset
centroids_k11 = kmeans_k11.cluster_centers_  # Get the centroids of the clusters
# Visualize the centroids for K=11 in a 1x11 grid and save the image as 'centroids_k11.png'
visualize_centroids(centroids_k11, 11, 'centroids_k11.png')

# This function calculates the clustering error by determining the majority label in each cluster. It counts the number of misclassified samples (those that don't match the majority) and sums them up across all clusters.
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
print(f"K=10 Error={error_k10}")

# Compute and print the clustering error for K=11
error_k11 = clustering_error(kmeans_k11, MNIST_y)
print(f"K=11 Error={error_k11}")
