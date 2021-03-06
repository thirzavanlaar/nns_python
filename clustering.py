#!/usr/bin/python

# apply hierarchical clustering to output cusize based on the 'Interaction Potential'

from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, centroid
import numpy as np
from scipy.spatial import distance
from haversine import haversine

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
cloudcentres = np.vstack((cloud_lon,cloud_lat)).T

# Compute distances for all pairs based on White et al 2018:
Y = distance.pdist(cloudcentres, haversine)

# Compute linkage matrix:
Z = centroid(Y)

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
plt.savefig('Figures/Dendrogram.pdf')

# Compute clusters, based on dendrogram and cut-off value max_d:
max_d = 1.0
clusters = fcluster(Z, max_d, criterion='distance')

# Plot clusters:
plt.figure(figsize=(10,8))
#plt.axis([-0.995, -0.991, 0.233, 0.236])
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.scatter(cloudcentres[:,0],cloudcentres[:,1], c=clusters)
plt.savefig('Figures/Clusters.pdf')


#--------------------------------------------------------------

rows = np.where(Z[:,3] == 2)
cloud1, cloud2, distance, count = Z[rows].T

histo = np.zeros(len(size))

for j in range(0,len(cloud1)):
    index1 = int(cloud1[j])
    index2 = int(cloud2[j])
    sizebin1 = int(cloud_bin[0,index1])
    sizebin2 = int(cloud_bin[0,index2])
    if sizebin1 < 4.:
        if sizebin2 >= 4.:
            histo[sizebin2] = histo[sizebin2] + 1
    if sizebin2 < 4.:
        if sizebin1 >= 4.:
            histo[sizebin1] = histo[sizebin1] + 1

print histo.shape

xvalues = np.arange(1, len(size) + 1)
print xvalues.shape

plt.figure()
plt.scatter(xvalues[0:25],histo[0:25])
plt.savefig('Figures/histogram.pdf')






