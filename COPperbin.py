#!/usr/bin/python

# Find minimal distances between clouds in one bin, average these per bin
# Compute geometric and arithmetical mean between all clouds per bin

from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
import numpy as np
from numpy import ma
from distance_methods import COP,COP_perbin

cusize = NetCDFFile(
         '/home/vanlaar/HDCP2data/TA_dom4/cusize_output_time41.nc')

size = cusize.variables['size']

begin_time = 41
end_time = 48


#######################################################

COP_all = np.zeros((end_time-begin_time+1,len(size)))

for time in range(begin_time,end_time+1):
    print 'time:',time

    cusize = NetCDFFile(
         '/home/vanlaar/HDCP2data/TA_dom4/cusize_output_time'+str(time)+'.nc')

    cloudlon = cusize.variables['cloud_lon']
    cloudlat = cusize.variables['cloud_lat']
    nclouds_real = cusize.variables['nclouds'][0]
    size = cusize.variables['size']
    cloudbin = cusize.variables['cloud_bin']
    ncloud_bin = cusize.variables['ncloud_bin'][:]

    ncloudsint = int(nclouds_real)
    cloud_lon = cloudlon[0,0:ncloudsint]
    cloud_lat = cloudlat[0,0:ncloudsint]
    cloud_bin = cloudbin[0,0:ncloudsint] 
    filledbin = int(max(cloud_bin))
    
    COP_output = COP_perbin(cloud_lon,cloud_lat,cloud_bin,size,filledbin)

    COP_all[time-41,:] = COP_output[:]


COP_perbin = np.mean(COP_all, axis=0)

####################################################################
#### Figures ####
####################################################################

JC = NetCDFFile('JosephCahalan.nc')
maxbin = JC.variables['maxbinvar'][:]


orange = (1.,0.38,0.01)
blue = (0.53,0.81,1)

plt.figure(figsize=(10,8))
plt.axis([0, 5500, 0, 0.15])
plt.xlabel('Cloud size', fontsize=15)
plt.ylabel('COP value', fontsize=15)
ax=plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.axvspan(size[maxbin], 5500, alpha =0.2, color='grey')
plt.scatter(size,COP_perbin)
plt.savefig('Figures/COP_perbin.pdf')

