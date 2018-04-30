#!/usr/bin/python

# apply hierarchical clustering to output cusize based on the 'Interaction Potential'

from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, centroid
import numpy as np
from scipy.spatial import distance
from haversine import haversine

cusize = NetCDFFile('cusize_output_subsub.nc')

# Import data from netcdf file:
cloudlon = cusize.variables['cloud_lon']
cloudlat = cusize.variables['cloud_lat']
nclouds_cusize  = cusize.variables['nclouds']
cloud_bin_netcdf = cusize.variables['cloud_bin']
size = cusize.variables['size']
nclouds = int(nclouds_cusize[0])
cloud_bin = np.array(cloud_bin_netcdf[0,0:nclouds:])


# Adjust data format for later use:
cloud_lon = cloudlon[0,0:nclouds]
cloud_lat = cloudlat[0,0:nclouds]
cloud_size = cloud_bin*size[0]
#cloudcentres = np.vstack((cloud_lon,cloud_lat,cloud_size)).T
cloudcentres = np.vstack((cloud_lon,cloud_lat)).T
labels = np.arange(nclouds)

# Compute distances for all pairs based on White et al 2018:
Y = distance.pdist(cloudcentres, haversine)

# Compute linkage matrix:
Z = centroid(Y)
max_d = 8000


# Plot dendrogram:
plt.figure(figsize=(25, 10))
plt.xlabel('Cloud label')
plt.ylabel('Euclidian distance [m]')
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
plt.savefig('dendrogram.pdf')

# Compute clusters, based on dendrogram and cut-off value max_d:
clusters = fcluster(Z, max_d, criterion='distance')

cloudcentres = np.rad2deg(cloudcentres)

# Plot clusters:
fig = plt.figure(figsize=(10,8))
plt.axis([-57, -56.8, 13.35, 13.5])
plt.xlabel('Longitude')
plt.ylabel('Latitude')
#ax = fig.add_subaxis(111)
#ax.imshow(world_map, extent[-57,-56.8,13.35,13.5], aspect = 'auto')
plt.scatter(cloudcentres[:,0],cloudcentres[:,1], c=clusters,s=cloud_size,cmap='prism')
for i, txt in enumerate(labels):
    plt.annotate(txt, (cloudcentres[i,0],cloudcentres[i,1]),
        xytext=(-5,-5),
        textcoords='offset points')

plt.savefig('clusters.pdf')


last = Z[-10:,2]
last_rev = last[::-1]
idxs = np.arange(1, len(last) + 1)
acceleration = np.diff(last,2)
acceleration_rev = acceleration[::-1]

plt.figure()
plt.plot(idxs,last_rev)
plt.plot(idxs[:-2] + 1, acceleration_rev)
plt.ylabel('Distance [m]')
plt.xlabel('Clusters')

plt.savefig('elbow.pdf')

#-------------------------------------------------------------

rows = np.where(Z[:,3] == 2)
cloud1, cloud2, distance, count = Z[rows].T

histo = np.zeros(len(size))

for j in range(0,len(cloud1)):
    index1 = int(cloud1[j])
    index2 = int(cloud2[j])
    sizebin1 = int(cloud_bin[index1])
    sizebin2 = int(cloud_bin[index2])
    if sizebin1 < 2.:
        histo[sizebin2] = histo[sizebin2] + 1
    if sizebin2 < 2.:
        histo[sizebin1] = histo[sizebin1] + 1
    
xvalues = np.arange(len(size))

plt.figure()
plt.scatter(xvalues[0:25],histo[0:25])
plt.savefig('histogram_subsub.pdf')

#--------------------------------------------------------------

print clusters

print max(cloud_size)





