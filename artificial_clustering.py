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
cloud_size2 = [2000,200,200,200,200,200,200,200,200]
labels = [0,1,2,3,4,5,6,7,8]
#cloudcentres = np.vstack((cloud_lon,cloud_lat,cloud_size)).T
cloudcentres = np.vstack((cloud_lon,cloud_lat)).T

# Compute distances for all pairs based on White et al 2018:
Y = distance.pdist(cloudcentres, euclidean)

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
max_d = 3
clusters = fcluster(Z, max_d, criterion='distance')

# Plot clusters:
plt.figure(figsize=(10,8))
plt.axis([-2, 10, -2, 10])
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.scatter(cloudcentres[:,0],cloudcentres[:,1], c=clusters,s=cloud_size2)
for i, txt in enumerate(labels):
    plt.annotate(txt, (cloudcentres[i,0],cloudcentres[i,1]),
        xytext=(-20,20),
        textcoords='offset points')
plt.savefig('Figures/Artificial_clusters.pdf')



#-----------------------------------------------------------------------
# Situation 2

size = [100,200]
nclouds = 7
labels = [0,1,2,3,4,5,6]
cloud_lon = [4,6,1,1,1,2,2]
cloud_lat = [2,6,5,6,7,6,7]
cloud_size = [2,2,1,1,1,1,1]
cloud_size2 = [2000,2000,200,200,200,200,200]
cloudcentres = np.vstack((cloud_lon,cloud_lat)).T
cloudcentres2 = np.vstack((cloud_lon,cloud_lat,cloud_size)).T

# Compute distances for all pairs based on White et al 2018:
Y = distance.pdist(cloudcentres, euclidean)

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
plt.axhline(1.2, c='black')
plt.savefig('Figures/Artificial_dendrogram_2A.pdf')

# Compute clusters, based on dendrogram and cut-off value max_d:
max_d = 1.2
clusters = fcluster(Z, max_d, criterion='distance')

# Plot clusters:
plt.figure(figsize=(10,8))
plt.axis([-2, 10, -2, 10])
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.scatter(cloudcentres[:,0],cloudcentres[:,1], c=clusters,s=cloud_size2)
for i, txt in enumerate(labels):
    plt.annotate(txt, (cloudcentres[i,0],cloudcentres[i,1]),
        xytext=(-20,20),
        textcoords='offset points')
plt.savefig('Figures/Artificial_clusters_2A.pdf')

# Now with V method
# Compute distances for all pairs based on White et al 2018:
Y = distance.pdist(cloudcentres2, euclidean_V)

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
plt.axhline(1,c='black')
plt.savefig('Figures/Artificial_dendrogram_2B.pdf')

# Compute clusters, based on dendrogram and cut-off value max_d:
max_d = 1
clusters = fcluster(Z, max_d, criterion='distance')

# Plot clusters:
plt.figure(figsize=(10,8))
plt.axis([-2, 10, -2, 10])
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.scatter(cloudcentres[:,0],cloudcentres[:,1], c=clusters,s=cloud_size2)
for i, txt in enumerate(labels):
    plt.annotate(txt, (cloudcentres[i,0],cloudcentres[i,1]),
        xytext=(-20,20),
        textcoords='offset points')
plt.savefig('Figures/Artificial_clusters_2B.pdf')

#-----------------------------------------------------------------------
# Situation 3

size = [100,200]
nclouds = 14
labels = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]
cloud_lon = [3,1,1,1,3,5,5,5,3,7,7,7,9,9,9]
cloud_lat = [4,2,4,6,6,6,4,2,2,4,2,6,2,4,6]
cloud_size = [2,1,1,1,1,1,1,1,1,2,1,1,1,1,1]
cloud_size2 = [2000,200,200,200,200,200,200,200,200,2000,200,200,200,200,200]
cloudcentres = np.vstack((cloud_lon,cloud_lat)).T
cloudcentres2 = np.vstack((cloud_lon,cloud_lat,cloud_size)).T

# Compute distances for all pairs based on White et al 2018:
Y = distance.pdist(cloudcentres, euclidean)

# Compute linkage matrix:
Z = ward(Y)

# Plot dendrogram:
plt.figure(figsize=(25, 10))
plt.xlabel('cloud label')
plt.ylabel('Euclidian distance')
dendrogram(
    Z,
    leaf_rotation=90.,
    #truncate_mode='lastp',
    #p=12,
    #show_leaf_counts=False,
    show_contracted=True,
)
plt.axhline(4, c='black')
plt.savefig('Figures/Artificial_dendrogram_3A.pdf')

# Compute clusters, based on dendrogram and cut-off value max_d:
max_d = 4
clusters = fcluster(Z, max_d, criterion='distance')

# Plot clusters:
plt.figure(figsize=(10,8))
plt.axis([-2, 10, -2, 10])
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.scatter(cloudcentres[:,0],cloudcentres[:,1], c=clusters,s=cloud_size2)
for i, txt in enumerate(labels):
    plt.annotate(txt, (cloudcentres[i,0],cloudcentres[i,1]),
        xytext=(-20,20),
        textcoords='offset points')
plt.savefig('Figures/Artificial_clusters_3A.pdf')

# Now with V method
# Compute distances for all pairs based on White et al 2018:
Y = distance.pdist(cloudcentres2, euclidean_V)

# Compute linkage matrix:
Z = ward(Y)

# Plot dendrogram:
plt.figure(figsize=(25, 10))
plt.xlabel('cloud label')
plt.ylabel('Euclidian distance')
dendrogram(
    Z,
    leaf_rotation=90.,
    #truncate_mode='lastp',
    #p=20,
    #show_leaf_counts=False,
    #show_contracted=True,
    #max_d = 0.8,
)
plt.axhline(0.8, c='black')
plt.savefig('Figures/Artificial_dendrogram_3B.pdf')

# Compute clusters, based on dendrogram and cut-off value max_d:
max_d = 0.8
clusters = fcluster(Z, max_d, criterion='distance')

# Plot clusters:
plt.figure(figsize=(10,8))
plt.axis([-2, 10, -2, 10])
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.scatter(cloudcentres[:,0],cloudcentres[:,1], c=clusters,s=cloud_size2)
for i, txt in enumerate(labels):
    plt.annotate(txt, (cloudcentres[i,0],cloudcentres[i,1]),
        xytext=(-20,20),
        textcoords='offset points')
plt.savefig('Figures/Artificial_clusters_3B.pdf')
