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
from distance_methods import distances

cusize = NetCDFFile(
         '/home/vanlaar/HDCP2data/TA_dom4/cusize_output_time41.nc')

size = cusize.variables['size']

begin_time = 41
end_time = 48

D0_all = np.zeros((end_time-begin_time+1,len(size)))
D1_all = np.zeros((end_time-begin_time+1,len(size)))
mindistance_mean_all = np.zeros((end_time-begin_time+1,len(size)))
mindistance_std_all = np.zeros((end_time-begin_time+1,len(size)))
hn_normalized_all = np.zeros((end_time-begin_time+1,len(size)))


for time in range(begin_time,end_time+1):
    print 'time:',time

    cusize = NetCDFFile(
         '/home/vanlaar/HDCP2data/TA_dom4/cusize_output_time'+str(time)+'.nc')

    cloudlon = cusize.variables['cloud_lon']
    cloudlat = cusize.variables['cloud_lat']
    nclouds_cusize  = cusize.variables['nclouds']
    #size = cusize.variables['size']
    cloud_bin = cusize.variables['cloud_bin']
    hn = cusize.variables['hn']
    hn_normalized_loop = hn/nclouds_cusize[0]

    nclouds = int(nclouds_cusize[0])
    cloud_lon = cloudlon[0,0:nclouds]
    cloud_lat = cloudlat[0,0:nclouds]
    filledbin=np.argmin(hn[0,:])  # last bin with clouds, rest is empty

    output_distances = distances(filledbin,cloud_lon,cloud_lat,cloud_bin,size,nclouds)

    D0_all[time-41] = output_distances[0]
    D1_all[time-41] = output_distances[1]
    mindistance_mean_all[time-41] = output_distances[2]
    mindistance_std_all[time-41] = output_distances[3]
    hn_normalized_all[time-41] = hn_normalized_loop

mindistance_mean = np.mean(mindistance_mean_all,axis=0)
mindistance_std = np.mean(mindistance_std_all,axis=0)
D0 = np.mean(D0_all,axis=0)
D1 = np.mean(D1_all,axis=0)
hn_normalized = np.mean(hn_normalized_all,axis=0)

sizelog = np.log10(size)
hnlog = ma.filled(np.log10(ma.masked_equal(hn_normalized, 0)), np.nan)

#res = ma.filled(log2(ma.masked_equal(m, 0)), 0)


mindistance_plus = mindistance_mean + mindistance_std
mindistance_minus = mindistance_mean - mindistance_std

filledbin = np.argmin(mindistance_mean)
slope, intercept, r_value, p_value, std_err = stats.linregress(size[0:filledbin],mindistance_mean[0:filledbin])

print 'r-squared:',r_value**2
line = intercept + slope*size

print 'slope:',slope
print 'intercept:',intercept


plt.figure(figsize=(10,8))
plt.axis([0, 5500, 0, 90000])
plt.xlabel('Cloud size [m]')
plt.ylabel('Nearest-neighbour distance [m]')
plt.fill_between(size,mindistance_plus,mindistance_minus,alpha=0.3,color='red')
plt.scatter(size,mindistance_mean,color='k')
#plt.scatter(size,mindistance_plus,color='g')
plt.plot(size,line,color='black')
plt.savefig('Figures/mindistance.pdf')


plt.figure(figsize=(10,8))
#plt.axis([50000,220000, 50000, 220000])
plt.xlabel('D1')
plt.ylabel('D0')
plt.scatter(D1,D0,color='k')
plt.savefig('Figures/D1-D0.pdf')


plt.figure(figsize=(10,8))
plt.xlabel('log(l) [m]')
plt.ylabel('log(N*(l)) [m-1]')
plt.scatter(sizelog,hnlog)
plt.savefig('Figures/CSD.pdf')


plt.figure()
plt.xlabel('size')
plt.ylabel('ratio distance/size')
plt.scatter(size[0:filledbin],mindistance_mean[0:filledbin]/size[0:filledbin])
plt.savefig('Figures/ratio_distance_size.pdf')
