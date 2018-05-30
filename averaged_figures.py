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
from CSD_fit import CSD_fit

cusize = NetCDFFile(
         '/home/vanlaar/HDCP2data/TA_dom4/cusize_output_time41.nc')

size = cusize.variables['size']

begin_time = 41
end_time = 48

D0_all = np.zeros((end_time-begin_time+1,len(size)))
D1_all = np.zeros((end_time-begin_time+1,len(size)))
nclouds_bin_all = np.zeros((end_time-begin_time+1,len(size)))
mindistance_mean_all = np.zeros((end_time-begin_time+1,len(size)))
mindistance_std_all = np.zeros((end_time-begin_time+1,len(size)))
maxdistance_all = np.zeros((end_time-begin_time+1,len(size)))
maxdistanceY_all = np.zeros((end_time-begin_time+1,len(size)))
hn_normalized_all = np.zeros((end_time-begin_time+1,len(size)))


for time in range(begin_time,end_time+1):
    print 'time:',time

    cusize = NetCDFFile(
         '/home/vanlaar/HDCP2data/TA_dom4/cusize_output_time'+str(time)+'.nc')

    cloudlon = cusize.variables['cloud_lon'][:]
    cloudlat = cusize.variables['cloud_lat'][:]
    nclouds_cusize  = cusize.variables['nclouds']
    #size = cusize.variables['size']
    cloud_bin = cusize.variables['cloud_bin'][0,:]
    hn = cusize.variables['hn']
    hn_normalized_loop = hn/nclouds_cusize[0]
    ncloud_bin = cusize.variables['ncloud_bin']

    ncloudsint = int(nclouds_cusize[0])
    cloud_lon = cloudlon[0,0:ncloudsint]
    cloud_lat = cloudlat[0,0:ncloudsint]
    filledbin=np.argmin(hn[0,:])  # last bin with clouds, rest is empty

    output_distances = distances(filledbin,cloud_lon,cloud_lat,cloud_bin,size,ncloudsint)

    D0_all[time-41] = output_distances[0]
    D1_all[time-41] = output_distances[1]
    mindistance_mean_all[time-41] = output_distances[2]
    mindistance_std_all[time-41] = output_distances[3]
    nclouds_bin_all[time-41] = output_distances[4]
    hn_normalized_all[time-41] = hn_normalized_loop


mindistance_mean = np.mean(mindistance_mean_all,axis=0)/1000
mindistance_std = np.mean(mindistance_std_all,axis=0)/1000
D0 = np.mean(D0_all,axis=0)
D1 = np.mean(D1_all,axis=0)
nclouds = np.mean(nclouds_bin_all,axis=0)
hn_normalized = np.mean(hn_normalized_all,axis=0)

filledbin_all=np.argmin(hn_normalized[:])  
fit = CSD_fit(hn_normalized[0:10],size[0:10])
logfit = fit[3]
a = fit[0]
b = fit[1]
c = fit[2]
print 'a, b, c:'
print a, b, c

print logfit

sizelog = np.log10(size)
hnlog = ma.filled(np.log10(ma.masked_equal(hn_normalized, 0)), np.nan)
ncloudslog = ma.filled(np.log10(ma.masked_equal(nclouds, 0)), np.nan)

#res = ma.filled(log2(ma.masked_equal(m, 0)), 0)

mindistance_plus = mindistance_mean + mindistance_std
mindistance_minus = mindistance_mean - mindistance_std

filledbin = np.argmin(mindistance_mean)
slope, intercept, r_value, p_value, std_err = stats.linregress(size[0:filledbin],mindistance_mean[0:filledbin])

print 'r-squared:',r_value**2
line = intercept + slope*size

print 'slope:',slope
print 'intercept:',intercept

##################################################################
### Plots

threshold = 0.005*sum(nclouds)
maxbin = np.min(np.where(nclouds <= threshold))

orange = (1.,0.38,0.01)
blue = (0.53,0.81,1)

plt.figure(figsize=(14,8))
plt.axis([0, 5500, 0, 120])
plt.xlabel('Cloud size [m]',fontsize=15)
plt.ylabel('Nearest-neighbour distance [km]',fontsize=15)
plt.fill_between(size,mindistance_plus,mindistance_minus,alpha=0.3,color=blue)
plt.scatter(size,mindistance_mean,color='k')
#plt.scatter(size,mindistance_plus,color='g')
#plt.plot(size,line,color='black')
ax = plt.gca()
ax.axvspan(size[maxbin], 5500, alpha=0.2, color='grey')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('Figures/mindistance.pdf')
plt.savefig('Figures/mindistance.png')


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
plt.scatter(sizelog[0:10],logfit[0:10])
plt.savefig('Figures/CSD.pdf')


plt.figure(figsize=(10,8))
plt.xlabel('Cloud size')
plt.ylabel('Ratio distance/size')
plt.axis([0, 5500, 0, 0.02])
ax = plt.gca()
ax.axvspan(size[maxbin], 5500, alpha=0.2, color='grey')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('Figures/mindistance.pdf')
plt.scatter(size[0:filledbin],mindistance_mean[0:filledbin]/size[0:filledbin])
plt.savefig('Figures/ratio_distance_size.pdf')



plt.figure(figsize=(10,8))
plt.xlabel('Cloud size')
plt.ylabel('Number of clouds')
plt.axhline(y=threshold, c='black')
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.scatter(size[0:filledbin],nclouds[0:filledbin])
plt.savefig('Figures/nclouds_size.pdf')

plt.figure(figsize=(10,8))
plt.xlabel('size')
plt.ylabel('nclouds')
plt.scatter(size[0:filledbin],ncloudslog[0:filledbin])
plt.savefig('Figures/nclouds_size_log.pdf')

