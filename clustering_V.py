#!/usr/bin/python

# apply hierarchical clustering to output cusize based on the 'Interaction Potential'

from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, ward
import numpy as np
from scipy.spatial import distance
from haversine import haversine_V

cusize = NetCDFFile('/home/vanlaar/HDCP2data/TA_dom4/cusize_output_time41.nc')

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
cloudcentres = np.vstack((cloud_lon,cloud_lat,cloud_size)).T

# Compute distances for all pairs based on White et al 2018:
Y = distance.pdist(cloudcentres, haversine_V)

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
plt.savefig('Figures/test_dendrogram.pdf')

# Compute clusters, based on dendrogram and cut-off value max_d:
max_d = 1.0
clusters = fcluster(Z, max_d, criterion='distance')

# Plot clusters:
plt.figure(figsize=(10,8))
#plt.axis([-0.995, -0.991, 0.233, 0.236])
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.scatter(cloudcentres[:,0],cloudcentres[:,1], c=clusters)
plt.savefig('Figures/test_clusters.pdf')









