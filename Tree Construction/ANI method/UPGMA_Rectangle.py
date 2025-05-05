import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt

# Step 1: Read the ANI similarity matrix and convert to a distance matrix
ani_df = pd.read_csv("./ani_output/ani_matrix.csv", index_col=0)
strains = ani_df.index.tolist()
ani_symmetric = (ani_df + ani_df.T) / 2

# Convert to distance matrix: distance = 1 - ANI
distance_matrix = 1.0 - ani_symmetric

# Verify symmetry (to prevent calculation errors)
assert np.allclose(distance_matrix, distance_matrix.T), "Distance matrix must be symmetric."

# SciPy requires the input to be in the form of compressed distance vectors (1D）
condensed_dist = squareform(distance_matrix)

# Step 2: Construct the hierarchical clustering tree
# Method options: "average" (UPGMA), "complete", "single", "ward", etc.
linkage_matrix = linkage(condensed_dist, method="average")

# Step 3: Visualizing tree structures (horizontal dendrograms)
plt.figure(figsize=(10, 5))
dendrogram(linkage_matrix, labels=strains, leaf_rotation=90)
plt.title("UPGMA(Average) Tree based on ANI distance (1 - ANI)")
plt.tight_layout()
plt.show()
