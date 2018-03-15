#!/usr/bin/python

# apply hierarchical clustering to output cusize based on the 'Interaction Potential'

from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, ward
import numpy as np
from scipy.spatial import distance
from artificial_distance import euclidean_V,euclidean


# Import data from netcdf file:
size = [100,200]
nclouds = 9

# Adjust data format for later use:
cloud_lon = [4,2,2,2,4,6,6,6,4]
cloud_lat = [4,2,4,6,6,6,4,2,2]
cloud_size = [2,1,1,1,1,1,1,1,1]
cloudcentres = np.vstack((cloud_lon,cloud_lat,cloud_size)).T

# Compute distances for all pairs based on White et al 2018:
Y = distance.pdist(cloudcentres, euclidean_V)

# Compute linkage matrix:
Z = ward(Y)

# Plot dendrogram:
plt.figure(figsize=(25, 10))
plt.xlabel('cloud label')
plt.ylabel('Euclidian distance')
dendrogram(
    Z,
    leaf_rotation=90.,
    truncate_mode='lastp',
    p=12,
    #show_leaf_counts=False,
    show_contracted=True,
)
plt.savefig('Figures/Artificial_dendrogram.pdf')

# Compute clusters, based on dendrogram and cut-off value max_d:
max_d = 1.4
clusters = fcluster(Z, max_d, criterion='distance')

# Plot clusters:
plt.figure(figsize=(10,8))
plt.axis([-2, 10, -2, 10])
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.scatter(cloudcentres[:,0],cloudcentres[:,1], c=clusters)
plt.savefig('Figures/Artificial_clusters.pdf')









