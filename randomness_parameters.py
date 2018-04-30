#!/usr/bin/python

# Find minimal distances between clouds in one bin, average these per bin
# Compute geometric and arithmetical mean between all clouds per bin

from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
import numpy as np
from numpy import ma
from distance_methods import distances,SCAI,COP
from Randomfield import randomfield

cusize = NetCDFFile(
         '/home/vanlaar/HDCP2data/TA_dom4/cusize_output_time41.nc')

size = cusize.variables['size']

grid = NetCDFFile('../NarvalDom2_NestE-R02B14_DOM03.nc')
clon_grid = grid.variables['clon']
clat_grid = grid.variables['clat']

Nmax = len(clon_grid)/2
lengthscale = 440000

begin_time = 41
end_time = 41

D0_all = np.zeros((end_time-begin_time+1,len(size)))
D1_all = np.zeros((end_time-begin_time+1,len(size)))
nclouds_bin_all = np.zeros((end_time-begin_time+1,len(size)))
nclouds_all = np.zeros((end_time-begin_time+1))
SCAI_0_all = np.zeros((end_time-begin_time+1))
SCAI_1_all = np.zeros((end_time-begin_time+1))

for time in range(begin_time,end_time+1):
    print 'time:',time

    cusize = NetCDFFile(
         '/home/vanlaar/HDCP2data/TA_dom4/cusize_output_time'+str(time)+'.nc')

    cloudlon = cusize.variables['cloud_lon']
    cloudlat = cusize.variables['cloud_lat']
    nclouds_cusize  = cusize.variables['nclouds']
    size = cusize.variables['size']
    cloudbin = cusize.variables['cloud_bin']
    hn = cusize.variables['hn']
    hn_normalized_loop = hn/nclouds_cusize[0]
    ncloud_bin = cusize.variables['ncloud_bin'][:]

    nclouds = int(nclouds_cusize[0])
    cloud_lon = cloudlon[0,0:nclouds]
    cloud_lat = cloudlat[0,0:nclouds]
    cloud_bin = cloudbin[0,0:nclouds] 
    filledbin = int(max(cloud_bin))
    
    output_distances = distances(filledbin,cloud_lon,cloud_lat,cloud_bin,size,nclouds)

    D0_all[time-41] = output_distances[0]
    D1_all[time-41] = output_distances[1]
    nclouds_bin_all[time-41] = output_distances[4]
    nclouds_all[time-41] = nclouds_cusize[0]

    #SCAIs = SCAI(cloud_lon,cloud_lat,nclouds,Nmax,lengthscale)

    #cloud_size = cloud_bin*size[0]
    #COP_output = COP(cloud_lon,cloud_lat,cloud_size)

    #SCAI_0_all[time-41] = SCAIs[0]
    #SCAI_1_all[time-41] = SCAIs[1]


D0 = np.mean(D0_all,axis=0)
D1 = np.mean(D1_all,axis=0)
nclouds_av = int(np.mean(nclouds_all))
nclouds_bin_av = np.mean(nclouds_bin_all,axis=0)
nrbins = len(nclouds_bin_av)

probability = nclouds_bin_av/nclouds_av

#print 'Real field:'
#print 'SCAI 0:', SCAI_0_all
#print 'SCAI 1:', SCAI_1_all
#print 'COP:', COP_output

#filledbin = np.argmin(mindistance_mean)
#slope, intercept, r_value, p_value, std_err = stats.linregress(size[0:filledbin],mindistance_mean[0:filledbin])
#line = intercept + slope*size

#############################################################
#### Randomness ####
############################################################

randomfield = randomfield(nclouds_av,nrbins,clon_grid,clat_grid,probability)

cloud_lon_random = randomfield[0]
cloud_lat_random = randomfield[1]
cloud_bin_random = np.array(randomfield[2])

#distances_random = distances(filledbin,cloud_lon_random,cloud_lat_random,cloud_bin_random,size,nrclouds)

#SCAIs = SCAI(cloud_lon_random,cloud_lat_random,nclouds_av,Nmax,lengthscale)
#
#SCAI_0_random = SCAIs[0]
#SCAI_1_random = SCAIs[1]
#
#cloud_size_random = cloud_bin_random*size[0]
#COP_random = COP(cloud_lon_random,cloud_lat_random,cloud_size_random)
#
#print 'Random:'
#print 'SCAI 0:', SCAI_0_random
#print 'SCAI 1:', SCAI_1_random
#print 'COP:', COP_random

unique, counts = np.unique(cloud_bin_random, return_counts=True)
frequency_random = np.array(counts).astype(np.float64)
frequency_randomT = frequency_random[:]/float(nclouds_av)

hn_random = np.zeros(nrbins,float)
hn_random[unique-1] = frequency_randomT[:]
print sum(hn_random)
hn_random_log = hn_random
hn_random_log[hn_random>0] = np.log10(hn_random[hn_random>0])

cloud_bin_real = nclouds_bin_av.astype(np.float64)
frequencyT = cloud_bin_real/float(nclouds_av)
print sum(frequencyT)
hn_real_log = frequencyT
hn_real_log[frequencyT>0] = np.log10(frequencyT[frequencyT>0])

size_log = np.log10(size)

plt.figure(figsize=(10,8))
plt.axis([2, 4.5, -4, -0.1])
plt.scatter(size_log,hn_real_log,c='k')
plt.scatter(size_log,hn_random_log,c='b')
plt.savefig('Figures/CSD_realrandom.pdf')


####################################################################
#### Figures ####
####################################################################

orange = (1.,0.38,0.01)
blue = (0.53,0.81,1)


plt.figure(figsize=(10,8))
#plt.axis([50000,220000, 50000, 220000])
plt.xlabel('D1')
plt.ylabel('D0')
plt.scatter(D1,D0,color='k')
plt.savefig('Figures/D1-D0.pdf')


plt.figure(figsize=(10,8))
plt.scatter(cloud_lon_random,cloud_lat_random,s=cloud_bin_random)
plt.savefig('randomfield.pdf')

