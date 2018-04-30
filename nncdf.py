#!/usr/bin/python

# Find minimal distances between clouds in one bin, average these per bin
# Compute geometric and arithmetical mean between all clouds per bin

from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
import numpy as np
from numpy import ma
from distance_methods import distances,NNCDF
from Randomfield import randomfield, randomfield_nooverlap

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
    ncloud_bin = cusize.variables['ncloud_bin']

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
    
    nncdf_real = NNCDF(cloud_lon,cloud_lat,size)


D0 = np.mean(D0_all,axis=0)
D1 = np.mean(D1_all,axis=0)
nclouds_av = int(np.mean(nclouds_all))
nclouds_bin_av = np.mean(nclouds_bin_all,axis=0)
nrbins = len(nclouds_bin_av)

probability = nclouds_bin_av/nclouds_av

#cloud_bin_av = np.zeros(nclouds_av)
#j=0 
#for i in range(0,nrbins):
#    elements = int(nclouds_bin_av[i])
#    if elements>0:
#        cloud_bin_av[j:j+elements] = i+1
#        j = j+elements
#        filledbin = i
#
#cloud_bin_av = np.trim_zeros(cloud_bin_av)
#nrclouds = len(cloud_bin_av)

#filledbin = np.argmin(mindistance_mean)
#slope, intercept, r_value, p_value, std_err = stats.linregress(size[0:filledbin],mindistance_mean[0:filledbin])
#line = intercept + slope*size

#############################################################
#### Randomness ####
############################################################

binwidth = size[0]

#randomfield = randomfield_nooverlap(nclouds_av,nrbins,clon_grid,clat_grid,probability,binwidth)
randomfield = randomfield(nclouds_av,nrbins,clon_grid,clat_grid,probability)

cloud_lon_random = randomfield[0]
cloud_lat_random = randomfield[1]
cloud_bin_random = np.array(randomfield[2])

#distances_random = distances(filledbin,cloud_lon_random,cloud_lat_random,cloud_bin_random,size,nrclouds)

nncdf_random = NNCDF(cloud_lon_random,cloud_lat_random,size)


plt.figure(figsize=(10,8))
plt.scatter(cloud_lon_random,cloud_lat_random,s=cloud_bin_random)
plt.savefig('Figures/randomfield.pdf')


plt.figure(figsize=(10,8))
plt.xlabel('nncdf random',fontsize=15)
plt.ylabel('nncdf real field',fontsize=15)
plt.scatter(nncdf_random,nncdf_real)
plt.plot(nncdf_random,nncdf_random,c='black')
plt.savefig('Figures/nncdf.pdf')

area = np.trapz(nncdf_real, nncdf_random)
print area



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


