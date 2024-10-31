# Dalya_Salih_Tier2.py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Set up colors and sizes according to the CPK coloring convention
atom_colors = {
    "C": "gray",   # Carbon
    "H": "white",  # Hydrogen
    "O": "red",    # Oxygen
    "N": "blue"    # Nitrogen
}
atom_sizes = {
    "C": 60,
    "H": 30,
    "O": 70,
    "N": 70
}

# 1. Load MDS coordinates and atom types
print("Loading MDS coordinates and atom types...")
molecule_coords = pd.read_csv("molecule_coordinates.csv")
molecule_distances = pd.read_csv("molecule_distances.tsv", sep='\t', header=0)

# Assume the atom types are in the first column
atom_types = molecule_distances.iloc[:, 0].values
print("Data loaded.")

# 2. Plotting in 3D with CPK coloring and bond connections
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.set_title('3D Molecular Structure')

# Iterate through each atom and plot it in 3D
for i, atom in molecule_coords.iterrows():
    atom_type = atom_types[i]
    color = atom_colors.get(atom_type, "black")  # default to black if atom type is unrecognized
    size = atom_sizes.get(atom_type, 50)         # default size if atom type is unrecognized
    ax.scatter(atom['X'], atom['Y'], atom['Z'], color=color, s=size, label=atom_type, alpha=0.7)

# 3. Draw bonds if the distance between atoms is < 1.6
for i in range(len(molecule_coords)):
    for j in range(i + 1, len(molecule_coords)):
        dist = np.linalg.norm(molecule_coords.iloc[i, :3] - molecule_coords.iloc[j, :3])
        if dist < 1.6:
            ax.plot([molecule_coords.iloc[i, 0], molecule_coords.iloc[j, 0]],
                    [molecule_coords.iloc[i, 1], molecule_coords.iloc[j, 1]],
                    [molecule_coords.iloc[i, 2], molecule_coords.iloc[j, 2]], color="black", linewidth=0.5)

# Save the plot
plt.savefig("molecule_3D_plot.png")
plt.show()
print("3D molecular plot saved as molecule_3D_plot.png.")
