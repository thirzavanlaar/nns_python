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
print len(time_input)
nr_timesteps = len(time_input)

deltat = datetime.timedelta(seconds=900)
print deltat
startdate = datetime.datetime(2013, 12, 20, 12, 05, 00, 00)
print startdate
#enddate = startdate + datetime.datetime(seconds=nr_timesteps*deltat)
enddate = startdate + datetime.timedelta(seconds=nr_timesteps*900)
print enddate

#base = datetime.datetime.today()
#date_list = [base - datetime.timedelta(days=x) for x in range(0, numdays)]

date_list = [startdate + datetime.timedelta(seconds=x*900) for x in range(0,nr_timesteps)]

print date_list[0]




parameters = NetCDFFile('random_parameters.nc')

SCAI_0_real = parameters.variables['SCAI_0_real'][:]
SCAI_1_real = parameters.variables['SCAI_1_real'][:]
COP_real = parameters.variables['COP_real'][:]

SCAI_0_random_min = parameters.variables['SCAI_0_random_min'][:]
SCAI_1_random_min = parameters.variables['SCAI_1_random_min'][:]
COP_random_min = parameters.variables['COP_random_min'][:]
Iorg_min = parameters.variables['Iorg_min'][:]

SCAI_0_random_max = parameters.variables['SCAI_0_random_max'][:]
SCAI_1_random_max = parameters.variables['SCAI_1_random_max'][:]
COP_random_max = parameters.variables['COP_random_max'][:]
Iorg_max = parameters.variables['Iorg_max'][:]


SCAI_0_relative_min = SCAI_0_real[:]/SCAI_0_random_min[:]
SCAI_1_relative_min = SCAI_1_real[:]/SCAI_1_random_min[:]
COP_relative_min = COP_real[:]/COP_random_min[:]
Iorg_relative_min = Iorg_min[:]/0.5

SCAI_0_relative_max = SCAI_0_real[:]/SCAI_0_random_max[:]
SCAI_1_relative_max = SCAI_1_real[:]/SCAI_1_random_max[:]
COP_relative_max = COP_real[:]/COP_random_max[:]
Iorg_relative_max = Iorg_max[:]/0.5

#xaxis = np.arange(0,len(COP_real))
xaxis = date_list[40:49]
print xaxis


#COP_adjusted = 1-(COP_relative-1)




plt.figure(figsize=(10,8))
plt.axis([date_list[39], enddate, 0.9, 1.25])
xfmt = md.DateFormatter('%H:%M:%S')
ax=plt.gca()
ax.xaxis.set_major_formatter(xfmt)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.xlabel('Time of day', fontsize=15)
plt.ylabel('Relative randomness', fontsize=15)
plt.scatter(xaxis,SCAI_0_relative_min, color='green', marker='o',label='SCAI 0, min')
plt.scatter(xaxis,SCAI_0_relative_max, color='green', marker='v',label='SCAI 0, max')
plt.scatter(xaxis,SCAI_1_relative_min, color='red', marker='o', label='SCAI 1, min')
plt.scatter(xaxis,SCAI_1_relative_max, color='red', marker='v', label='SCAI 1, max')
plt.scatter(xaxis,COP_relative_min, color='black', marker='o', label='COP, min')
plt.scatter(xaxis,COP_relative_max, color='black', marker='v', label='COP, max')
plt.scatter(xaxis,Iorg_relative_min, color='blue', marker='o', label='Iorg, min')
plt.scatter(xaxis,Iorg_relative_max, color='blue', marker='v', label='Iorg, max')
plt.axhline(1, c='k')
plt.legend()
plt.savefig('Figures/random_parameter_minmax.pdf')


