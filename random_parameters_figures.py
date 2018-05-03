#!/usr/bin/python

# Find minimal distances between clouds in one bin, average these per bin
# Compute geometric and arithmetical mean between all clouds per bin

from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
import numpy as np


parameters = NetCDFFile('random_parameters.nc')

SCAI_0_real = parameters.variables['SCAI_0_real'][:]
SCAI_1_real = parameters.variables['SCAI_1_real'][:]
COP_real = parameters.variables['COP_real'][:]

SCAI_0_random = parameters.variables['SCAI_0_random'][:]
SCAI_1_random = parameters.variables['SCAI_1_random'][:]
COP_random = parameters.variables['COP_random'][:]

Iorg = parameters.variables['Iorg'][:]

SCAI_0_relative = SCAI_0_real[:]/SCAI_0_random[:]
SCAI_1_relative = SCAI_1_real[:]/SCAI_1_random[:]
COP_relative = COP_real[:]/COP_random[:]
Iorg_relative = Iorg[:]/0.5

xaxis = np.arange(0,len(COP_real))

COP_adjusted = 1-(COP_relative-1)


plt.figure(figsize=(10,8))
plt.xlabel('Snapshot nr', fontsize=15)
plt.ylabel('Relative randomness', fontsize=15)
plt.scatter(xaxis,SCAI_0_relative, label='SCAI 0')
plt.scatter(xaxis,SCAI_1_relative, label='SCAI 1')
plt.scatter(xaxis,COP_adjusted, label='COP')
plt.scatter(xaxis,Iorg_relative, label='Iorg')
plt.axhline(1, c='k')
plt.legend()
plt.savefig('Figures/random_parameter.pdf')


