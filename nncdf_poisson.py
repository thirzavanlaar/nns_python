#!/usr/bin/python

# Find minimal distances between clouds in one bin, average these per bin
# Compute geometric and arithmetical mean between all clouds per bin

from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
import numpy as np
import time
import datetime
from matplotlib.dates import drange
import matplotlib.dates as md
from math import exp, factorial, pi

cusize = NetCDFFile('/home/vanlaar/HDCP2data/TA_dom4/python_output/cusize_output.nc')

size = cusize.variables['size'][:]

parameters = NetCDFFile('random_averaged_test.nc')

SCAI_0_real = parameters.variables['SCAI_0_real'][:]
SCAI_1_real = parameters.variables['SCAI_1_real'][:]
COP_real = parameters.variables['COP_real'][:]

nclouds = parameters.variables['nclouds'][40]
nncdf_real = parameters.variables['nncdf'][40,:]
bins = parameters.variables['nncdf_binedges'][40,:]

print nncdf_real*100


domainsize = 6e10
density = nclouds/domainsize

#poisson = (nclouds**binnr[40]*exp(-1*nclouds))/factorial(binnr[40])

r = bins

poisson = np.zeros(len(r))

nr = 0
for element in r:
    poisson[nr] = 1-exp(-density*pi*element**2)
    nr += 1


xaxis = np.arange(0,99)
xaxis2 = np.arange(0,100)


plt.figure(figsize=(10,8))
#plt.axis([date_list[39], enddate, 0.9, 1.25])
#xfmt = md.DateFormatter('%H:%M:%S')
ax=plt.gca()
#ax.xaxis.set_major_formatter(xfmt)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.xlabel('distance', fontsize=15)
plt.ylabel('nncdf', fontsize=15)
plt.scatter(xaxis,nncdf_real*100, color='green', marker='o', label='real')
plt.scatter(xaxis2,poisson, color='red', marker='o', label='poisson')
#plt.axhline(1, c='k')
plt.legend()
plt.savefig('Figures/nncdf_poisson.pdf')


