#!/usr/bin/python

# Find minimal distances between clouds in one bin, average these per bin
# Compute geometric and arithmetical mean between all clouds per bin

from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats 
from haversine import haversine
from scipy.spatial import distance

# Import data:

cusize = NetCDFFile(
         '/home/vanlaar/HDCP2data/TA_dom4/cusize_output_time41.nc')

cloudlon = cusize.variables['cloud_lon']
cloudlat = cusize.variables['cloud_lat']
nclouds_cusize  = cusize.variables['nclouds']
size = cusize.variables['size']
cloud_bin = cusize.variables['cloud_bin']
hn = cusize.variables['hn']

nclouds = int(nclouds_cusize[0])
cloud_lon = cloudlon[0,0:nclouds]
cloud_lat = cloudlat[0,0:nclouds]
filledbin=np.argmin(hn[0,:])  # last bin with clouds, rest is empty



def distances(filledbin,cloud_lon,cloud_lat,cloud_bin):
    ncloud_bin=np.zeros((len(size)))
    D1=np.zeros((len(size)))
    D0=np.zeros((len(size)))
    mindistance_mean=np.zeros((len(size)))
    mindistance_std=np.zeros((len(size)))

    for bb in range(0, filledbin):
        binclouds=np.zeros((nclouds,2))
        # Select clouds present in current bin bb:
        idx = np.where(cloud_bin[0,:]==bb+1)          
        binclouds[idx,0]=cloud_lon[idx]
        binclouds[idx,1]=cloud_lat[idx]
        ncloud_bin[bb] = len(idx[0])
        # Compute all distances by using the haversine method:
        Y = distance.pdist(binclouds[idx],haversine)
        D0[bb] = stats.gmean(Y)
        D1[bb] = np.mean(Y)

        # Taking only the minimum distance per cloud (nearest neighbours):
        Z = distance.squareform(Y)
        Z[Z==0] = np.nan
        mindistances = np.nanmin(Z,axis=0)
        mindistance_mean[bb] = np.mean(mindistances)
        mindistance_std[bb] = np.std(mindistances)
    return D0,D1,mindistance_mean,mindistance_std,ncloud_bin
        

output_distances = distances(filledbin,cloud_lon,cloud_lat,cloud_bin)
D0 = output_distances[0]
D1 = output_distances[1]
mindistance_mean = output_distances[2]
mindistance_std = output_distances[3]


mindistance_plus = mindistance_mean+mindistance_std

plt.figure(figsize=(10,8))
plt.axis([0, 6000, 0, 70000])
plt.xlabel('Size')
plt.ylabel('Minimal distance')
plt.scatter(size,mindistance_mean,color='k')
plt.scatter(size,mindistance_plus,color='g')
plt.savefig('Figures/mindistance.pdf')


plt.figure(figsize=(10,8))
#plt.axis([50000,220000, 50000, 220000])
plt.xlabel('D1')
plt.ylabel('D0')
plt.scatter(D1,D0,color='k')
plt.savefig('Figures/D1-D0.pdf')

