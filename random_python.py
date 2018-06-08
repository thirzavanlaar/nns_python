#!/usr/bin/python

# Find minimal distances between clouds in one bin, average these per bin
# Compute geometric and arithmetical mean between all clouds per bin

from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
import numpy as np
from Randomfield import randomfield, randomfield_nooverlap, overlap

# Grid info:

grid = NetCDFFile('NarvalDom2_NestE-R02B14_DOM03.nc')
clon_grid = grid.variables['clon']
clat_grid = grid.variables['clat']

# cusize variables:

cusize = NetCDFFile('/home/vanlaar/HDCP2data/TA_dom4/python_output/cusize_output_old.nc')

nclouds = cusize.variables['nclouds'][40]
ncloudsint = int(nclouds)
cloudlon = cusize.variables['cloud_lon'][40,0:ncloudsint]
cloudlat = cusize.variables['cloud_lat'][40,0:ncloudsint]
cloudsize = cusize.variables['cloud_size'][40,0:ncloudsint]
size = cusize.variables['size'][:]
CSD = cusize.variables['CSD'][40,:]

print nclouds

# for random field:

#print CSD
#frequency = CSD*nclouds
#print frequency
#probability = frequency/nclouds
#print probability
#print sum(CSD)

nrbins = len(size)

minbinsize = size[0]
cloud_bin = np.ceil(cloudsize[:]/minbinsize)

unique, counts = np.unique(cloud_bin, return_counts=True)
indeces = np.array(unique).astype(np.int)

frequency = np.zeros(nrbins,float)
frequency[indeces-1] = counts

probability = frequency/nclouds





randomfield = randomfield(ncloudsint,nrbins,clon_grid,clat_grid,probability)

cloud_lon_random = randomfield[0]
cloud_lat_random = randomfield[1]
cloud_bin_random = np.array(randomfield[2])


plt.figure(figsize=(10,8))
#plt.axis([-0.995, -0.991, 0.233, 0.236])
plt.scatter(cloud_lon_random,cloud_lat_random,s=cloud_bin_random)
plt.savefig('Figures/random_python.pdf')

# now without overlap

binwidth = size[0]

randomfield_no = randomfield_nooverlap(ncloudsint,nrbins,clon_grid,clat_grid,probability,binwidth)

cloud_lon_random_no = randomfield_no[0]
cloud_lat_random_no = randomfield_no[1]
cloud_bin_random_no = np.array(randomfield_no[2])



#randomfield_overlap = overlap(ncloudsint,nrbins,clon_grid,clat_grid,probability,binwidth)
#
#cloud_lon_random_no = randomfield_overlap[0]
#cloud_lat_random_no = randomfield_overlap[1]
#cloud_bin_random_no = np.array(randomfield_overlap[2])
#overlap_no = randomfield_overlap[3]
#
#print overlap_no

