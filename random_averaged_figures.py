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


time_info = NetCDFFile('/data/inscape/phil/HDCP2/HDCP2TA_RegE_PR1250m-P01_CLD_DOM04_ML_20131220T120000Z_qc.nc')

time_input = time_info.variables['time'][:]
nr_timesteps = len(time_input)

deltat = datetime.timedelta(seconds=900)
startdate = datetime.datetime(2013, 12, 20, 12, 05, 00, 00)
#enddate = startdate + datetime.datetime(seconds=nr_timesteps*deltat)
enddate = startdate + datetime.timedelta(seconds=nr_timesteps*900)

#base = datetime.datetime.today()
#date_list = [base - datetime.timedelta(days=x) for x in range(0, numdays)]

date_list = [startdate + datetime.timedelta(seconds=x*900) for x in range(0,nr_timesteps)]

print date_list[0]


parameters = NetCDFFile('random_averaged.nc')

SCAI_0_real = parameters.variables['SCAI_0_real'][:]
SCAI_1_real = parameters.variables['SCAI_1_real'][:]
COP_real = parameters.variables['COP_real'][:]
nclouds = parameters.variables['nclouds'][:]


#xaxis = np.arange(0,len(COP_real))
xaxis = date_list[0:48]
#print xaxis


#COP_adjusted = 1-(COP_relative-1)

print len(xaxis)
print SCAI_0_real.shape


plt.figure(figsize=(10,8))
plt.axis([date_list[0], enddate, 0.0, 1.4])
xfmt = md.DateFormatter('%H:%M:%S')
ax=plt.gca()
ax.xaxis.set_major_formatter(xfmt)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.xlabel('Time of day', fontsize=15)
plt.ylabel('Randomness', fontsize=15)
plt.scatter(xaxis,SCAI_0_real, color='green', marker='o',label='SCAI 0')
plt.scatter(xaxis,SCAI_1_real, color='red', marker='o', label='SCAI 1')
#plt.axhline(1, c='k')
plt.legend()
plt.savefig('Figures/SCAI_timeevolution.pdf')

plt.figure(figsize=(10,8))
plt.axis([date_list[0], enddate, 0.01, 0.03])
xfmt = md.DateFormatter('%H:%M:%S')
ax=plt.gca()
ax.xaxis.set_major_formatter(xfmt)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.xlabel('Time of day', fontsize=15)
plt.ylabel('Randomness', fontsize=15)
plt.scatter(xaxis,COP_real, color='blue', marker='o', label='COP')
#plt.axhline(1, c='k')
plt.legend()
plt.savefig('Figures/COP_timeevolution.pdf')

plt.figure(figsize=(10,8))
plt.axis([date_list[0], enddate, 0, 6000])
xfmt = md.DateFormatter('%H:%M:%S')
ax=plt.gca()
ax.xaxis.set_major_formatter(xfmt)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.xlabel('Time of day', fontsize=15)
plt.ylabel('nclouds', fontsize=15)
plt.scatter(xaxis,nclouds, color='black', marker='o',label='nclouds')
#plt.axhline(1, c='k')
plt.legend()
plt.savefig('Figures/nclouds_timeevolution.pdf')



