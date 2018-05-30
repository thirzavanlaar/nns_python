#!/usr/bin/python

# Find minimal distances between clouds in one bin, average these per bin
# Compute geometric and arithmetical mean between all clouds per bin

from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
import numpy as np
from numpy import ma

JC = NetCDFFile('JosephCahalan.nc')

size = JC.variables['size'][:]
mindistance_min = JC.variables['mindistance_min'][:]
mindistance_plus = JC.variables['mindistance_plus'][:]
mindistance_mean = JC.variables['mindistance_mean'][:]
neighbour_size = JC.variables['neighbour_size'][:]
maxbin = JC.variables['maxbinvar'][:]
filledbin = JC.variables['filledbinvar'][:]

print maxbin



begin_time = 41
end_time = 48


orange = (0.93,0.47,0.26)
blue = (0.53,0.81,1)
green = (0.13,0.55,0.13)


x_JC = size[0:5]
y_JC = -0.38 + 16*x_JC 





plt.figure(figsize=(14,8))
plt.axis([0, 5500, 0, 4])
plt.xlabel('Cloud size [m]',fontsize=15)
plt.ylabel('Nearest-neighbour distance [km]',fontsize=15)
plt.fill_between(size,mindistance_plus,mindistance_min,alpha=0.3,color=blue)
plt.scatter(size,mindistance_mean,color='k')
#plt.plot(x_JC, y_JC/1000)
#plt.scatter(size,mindistance_plus,color='g')
#plt.plot(size,line,color='black')
ax = plt.gca()
ax.axvspan(size[maxbin], 5500, alpha=0.2, color='grey')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('Figures/mindistance_JC.pdf')
plt.savefig('Figures/mindistance_JC.png')



plt.figure(figsize=(14,8))
plt.axis([0, 5500, 0, 0.008])
plt.xlabel('Cloud size',fontsize=15)
plt.ylabel('Ratio distance/size',fontsize=15)
plt.scatter(size[0:filledbin],mindistance_mean[0:filledbin]/size[0:filledbin])
ax = plt.gca()
ax.axvspan(size[maxbin], 5500, alpha=0.2, color='grey')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('Figures/ratio_distance_size_JC.pdf')


plt.figure(figsize=(14,8))
plt.xlabel('Cloud size [m]',fontsize=15)
plt.ylabel('Nearest neighbour size [m]',fontsize=15)
plt.axis([0, 5500, 0, 800])
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.axvspan(size[maxbin], 5500, alpha=0.2, color='grey')
plt.scatter(size[0:filledbin],neighbour_size[0:filledbin],c=green,s=120)
plt.savefig('Figures/neighbour_size_avg.pdf')
plt.savefig('Figures/neighbour_size_avg.png')

