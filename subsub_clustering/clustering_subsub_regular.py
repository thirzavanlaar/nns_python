#!/usr/bin/python

# apply hierarchical clustering to output cusize based on the 'Interaction Potential'

from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, ward
import numpy as np
from scipy.spatial import distance
from haversine import haversine

cusize = NetCDFFile('cusize_output_subsub_regular.nc')

# Import data from netcdf file:
cloudlon = cusize.variables['cloud_lon']
cloudlat = cusize.variables['cloud_lat']
nclouds_cusize  = cusize.variables['nclouds']
cloud_bin = cusize.variables['cloud_bin']
size = cusize.variables['size']
nclouds = int(nclouds_cusize[0])

# Adjust data format for later use:
cloud_lon = cloudlon[0,0:nclouds]
cloud_lat = cloudlat[0,0:nclouds]
cloud_size = cloud_bin[0,:nclouds]*size[0]
#cloudcentres = np.vstack((cloud_lon,cloud_lat,cloud_size)).T
cloudcentres = np.vstack((cloud_lon,cloud_lat)).T

# Compute distances for all pairs based on White et al 2018:
Y = distance.pdist(cloudcentres, haversine)

# Compute linkage matrix:
Z = ward(Y)
max_d = 4000

# Plot dendrogram:
plt.figure(figsize=(25, 10))
plt.xlabel('cloud label')
plt.ylabel('Euclidian distance')
dendrogram(
    Z,
    leaf_rotation=90.,
    truncate_mode='lastp',
    p=50,
    #show_leaf_counts=False,
    show_contracted=True,
    color_threshold=max_d,
)
plt.axhline(max_d, c='black')
plt.savefig('dendrogram_regular.pdf')

# Compute clusters, based on dendrogram and cut-off value max_d:
clusters = fcluster(Z, max_d, criterion='distance')

cloudcentres = np.rad2deg(cloudcentres)

# Plot clusters:
plt.figure(figsize=(10,8))
#plt.axis([-0.995, -0.991, 0.233, 0.236])
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.scatter(cloudcentres[:,0],cloudcentres[:,1], c=clusters,s=cloud_size,cmap='prism')
plt.savefig('clusters_regular.pdf')









