#!/usr/bin/python

# apply hierarchical clustering to output cusize

from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import numpy as np

cusize = NetCDFFile('/home/vanlaar/HDCP2data/TA_dom4/cusize_output_time41.nc')

cloudlon = cusize.variables['cloud_lon']
cloudlat = cusize.variables['cloud_lat']
nclouds_cusize  = cusize.variables['nclouds']

nclouds = int(nclouds_cusize[0])

cloud_lon = cloudlon[0,0:nclouds]
cloud_lat = cloudlat[0,0:nclouds]


#length = len(cloud_lon[0,:])

cloudcentres = np.column_stack((cloud_lon,cloud_lat))

Z = linkage(cloudcentres,'ward')


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
plt.savefig('dendrogram.pdf')
#plt.close(fig)

max_d = 0.1
clusters = fcluster(Z, max_d, criterion='distance')

plt.figure(figsize=(10,8))
#plt.axis([-0.995, -0.991, 0.233, 0.236])
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.scatter(cloudcentres[:,0],cloudcentres[:,1], c=clusters)
plt.savefig('clusters.pdf')
#plt.close(fig)









