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
from distance_methods import distances,JosephCahalan

cusize = NetCDFFile(
         '/home/vanlaar/HDCP2data/TA_dom4/cusize_output_time41.nc')

size = cusize.variables['size']

begin_time = 41
end_time = 48

mindistance_JC_mean_all = np.zeros((end_time-begin_time+1,len(size)))
mindistance_JC_std_all = np.zeros((end_time-begin_time+1,len(size)))


for time in range(begin_time,end_time+1):
    print 'time:',time

    cusize = NetCDFFile(
         '/home/vanlaar/HDCP2data/TA_dom4/cusize_output_time'+str(time)+'.nc')

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

    output = JosephCahalan(filledbin,cloud_bin,cloud_lon,cloud_lat,size)

    mindistance_JC_mean_all[time-41] = output[0]
    mindistance_JC_std_all[time-41] = output[1]

mindistance_JC_mean = np.mean(mindistance_JC_mean_all,axis=0)
mindistance_JC_std = np.mean(mindistance_JC_std_all,axis=0)

mindistance_JC_plus = mindistance_JC_mean + mindistance_JC_std
mindistance_JC_minus = mindistance_JC_mean - mindistance_JC_std

filledbin = np.argmin(mindistance_JC_mean)
slope, intercept, r_value, p_value, std_err = stats.linregress(size[0:filledbin],mindistance_JC_mean[0:filledbin])

print 'r-squared:',r_value**2
line = intercept + slope*size

print 'slope:',slope
print 'intercept:',intercept


plt.figure(figsize=(10,8))
plt.axis([0, 5500, 0, 4000])
plt.xlabel('Cloud size [m]')
plt.ylabel('Nearest-neighbour distance [m]')
plt.fill_between(size,mindistance_JC_plus,mindistance_JC_minus,alpha=0.3,color='red')
plt.scatter(size,mindistance_JC_mean,color='k')
#plt.scatter(size,mindistance_plus,color='g')
plt.plot(size,line,color='black')
plt.savefig('Figures/mindistance_JC.pdf')



plt.figure()
plt.xlabel('size')
plt.ylabel('ratio distance/size')
plt.scatter(size[0:filledbin],mindistance_JC_mean[0:filledbin]/size[0:filledbin])
plt.savefig('Figures/ratio_distance_size_JC.pdf')
