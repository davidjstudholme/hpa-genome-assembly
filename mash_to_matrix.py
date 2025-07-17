import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Input/output files from command-line args
dist_file = sys.argv[1]
matrix_out = sys.argv[2]
heatmap_out = sys.argv[3]

# Read Mash output
df = pd.read_csv(dist_file, sep="\t", header=None,
                 names=["genome1", "genome2", "distance", "pvalue", "matches"])

# Build a symmetric matrix
matrix = df.pivot(index="genome1", columns="genome2", values="distance")
matrix = matrix.combine_first(matrix.T)
matrix = matrix.fillna(0)

# Save as CSV
matrix.to_csv(matrix_out)

# Plot heatmap
sns.set(font_scale=0.8)
sns.clustermap(matrix, cmap="viridis", annot=True, fmt=".3f", figsize=(10, 10))
plt.savefig(heatmap_out, dpi=300)
