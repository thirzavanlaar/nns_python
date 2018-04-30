#!/usr/bin/python

# Find minimal distances between clouds in one bin, average these per bin
# Compute geometric and arithmetical mean between all clouds per bin

from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
import numpy as np
from numpy import ma
from scipy import stats 
from haversine import haversine
from scipy.spatial import distance
from distance_methods import distances, JosephCahalan
import matplotlib as mpl

cusize = NetCDFFile(
         '/home/vanlaar/HDCP2data/TA_dom4/cusize_output_time41.nc')


size = cusize.variables['size']
cloudlon = cusize.variables['cloud_lon']
cloudlat = cusize.variables['cloud_lat']
nclouds_cusize  = cusize.variables['nclouds']
cloud_bin = cusize.variables['cloud_bin']  # binnr of clouds
hn = cusize.variables['hn']
hn_normalized = hn/nclouds_cusize[0]

nclouds = int(nclouds_cusize[0])
cloud_lon = cloudlon[0,0:nclouds]
cloud_lat = cloudlat[0,0:nclouds]
filledbin=np.argmin(hn[0,:])  # last bin with clouds, rest is empty

output_distances = JosephCahalan(filledbin,cloud_bin,cloud_lon,cloud_lat,size)

mindistance_JC_mean = output_distances[0]
mindistance_JC_std = output_distances[1]
neighbour_avg = output_distances[2]
neighbour_histo = output_distances[3]


sizelog = np.log10(size)
#hn_normalized = np.ma.masked_where(hn_normalized==0.,hn_normalized)
#hnlog = np.log10(hn_normalized)

#plt.figure(figsize=(10,8))
#plt.hist(neighbour_size)
#plt.show()

end = 30
xaxis = np.arange(len(size))

neighbour_histo = np.ma.masked_where(neighbour_histo==0,neighbour_histo)

cmap = mpl.cm.gnuplot


plt.figure(figsize=(10,8))
for i in range(0,filledbin):
    plt.plot(xaxis[0:end],neighbour_histo[i,0:end],color=cmap(i/float(filledbin)))
plt.savefig('Figures/histogram_TA.pdf')

print neighbour_avg.shape

plt.figure(figsize=(10,8))
plt.scatter(xaxis[0:filledbin],neighbour_avg[0:filledbin])
plt.savefig('Figures/neighbour_size.pdf')

