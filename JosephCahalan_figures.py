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
from distance_methods import distances,JosephCahalan, JosephCahalan_kneighbour

cusize = NetCDFFile(
         '/home/vanlaar/HDCP2data/TA_dom4/cusize_output_time41.nc')

size = cusize.variables['size']

begin_time = 41
end_time = 48

mindistance_JC_mean_all = np.zeros((end_time-begin_time+1,len(size)))
mindistance_JC_std_all = np.zeros((end_time-begin_time+1,len(size)))
neighbour_avg_all = np.zeros((end_time-begin_time+1,len(size)))


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

    #output = JosephCahalan(filledbin,cloud_bin,cloud_lon,cloud_lat,size)
    output = JosephCahalan_kneighbour(filledbin,cloud_bin,cloud_lon,cloud_lat,size,5)

    mindistance_JC_mean_all[time-41] = output[0]
    mindistance_JC_std_all[time-41] = output[1]
    neighbour_avg_all[time-41,:] = output[2]


mindistance_JC_mean = np.mean(mindistance_JC_mean_all,axis=0)/1000
mindistance_JC_std = np.mean(mindistance_JC_std_all,axis=0)/1000
neighbour_avg = np.mean(neighbour_avg_all,axis=0)*size[0]

mindistance_JC_plus = mindistance_JC_mean + mindistance_JC_std
mindistance_JC_minus = mindistance_JC_mean - mindistance_JC_std

filledbin = np.argmin(mindistance_JC_mean)
slope, intercept, r_value, p_value, std_err = stats.linregress(size[0:filledbin],mindistance_JC_mean[0:filledbin])

print 'r-squared:',r_value**2
line = intercept + slope*size

print 'slope:',slope
print 'intercept:',intercept

orange = (0.93,0.47,0.26)
blue = (0.53,0.81,1)
green = (0.13,0.55,0.13)

plt.figure(figsize=(14,8))
plt.axis([0, 5500, 0, 4])
plt.xlabel('Cloud size [m]',fontsize=15)
plt.ylabel('Nearest-neighbour distance [km]',fontsize=15)
plt.fill_between(size,mindistance_JC_plus,mindistance_JC_minus,alpha=0.3,color=blue)
plt.scatter(size,mindistance_JC_mean,color='k')
#plt.scatter(size,mindistance_plus,color='g')
#plt.plot(size,line,color='black')
plt.savefig('Figures/mindistance_JC.pdf')
plt.savefig('Figures/mindistance_JC.png')



plt.figure()
plt.xlabel('size',fontsize=15)
plt.ylabel('ratio distance/size',fontsize=15)
plt.scatter(size[0:filledbin],mindistance_JC_mean[0:filledbin]/size[0:filledbin])
plt.savefig('Figures/ratio_distance_size_JC.pdf')


plt.figure(figsize=(14,8))
plt.xlabel('Cloud size [m]',fontsize=15)
plt.ylabel('Nearest neighbour size [m]',fontsize=15)
plt.scatter(size[0:filledbin],neighbour_avg[0:filledbin],c=green,s=120)
plt.savefig('Figures/neighbour_size_avg.pdf')
plt.savefig('Figures/neighbour_size_avg.png')

