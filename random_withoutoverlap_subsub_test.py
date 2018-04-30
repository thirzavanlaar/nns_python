#!/usr/bin/python

# Find minimal distances between clouds in one bin, average these per bin
# Compute geometric and arithmetical mean between all clouds per bin

from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
import numpy as np
from Randomfield import randomfield, randomfield_nooverlap

# Grid info:

grid = NetCDFFile('subsub_clustering/subsubgrid.nc')
clon_grid = grid.variables['clon']
clat_grid = grid.variables['clat']

# cusize variables:

cusize = NetCDFFile('cusize_output_subsub.nc')

cloudlon = cusize.variables['cloud_lon']
cloudlat = cusize.variables['cloud_lat']
nclouds = cusize.variables['nclouds'][:]
cloud_bin_netcdf = cusize.variables['cloud_bin']
size = cusize.variables['size']
ncloudsint = int(nclouds[0])
cloud_bin = np.array(cloud_bin_netcdf[0,0:ncloudsint:])

cloud_lon = cloudlon[0,0:ncloudsint]
cloud_lat = cloudlat[0,0:ncloudsint]
cloud_size = cloud_bin*size[0]

# for random field:

unique, counts = np.unique(cloud_bin, return_counts=True)
nrbins = len(size)
indeces = np.array(unique).astype(np.int)

frequency = np.zeros(nrbins,float)
frequency[indeces-1] = counts

probability = frequency/nclouds

print nclouds

ncloudsint = 50

randomfield = randomfield(ncloudsint,nrbins,clon_grid,clat_grid,probability)

cloud_lon_random = randomfield[0]
cloud_lat_random = randomfield[1]
cloud_bin_random = np.array(randomfield[2])


plt.figure(figsize=(10,8))
plt.axis([-0.995, -0.991, 0.233, 0.236])
plt.scatter(cloud_lon_random,cloud_lat_random,s=cloud_bin_random)
plt.savefig('Figures/random_subsub.pdf')

# no without overlap

binwidth = size[0]

randomfield_no = randomfield_nooverlap(ncloudsint,nrbins,clon_grid,clat_grid,probability,binwidth)

cloud_lon_random_no = randomfield_no[0]
cloud_lat_random_no = randomfield_no[1]
cloud_bin_random_no = np.array(randomfield_no[2])




