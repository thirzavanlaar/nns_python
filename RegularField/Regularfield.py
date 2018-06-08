#!/usr/bin/python


from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
import numpy as np

grid = NetCDFFile('subsubgrid.nc')

clon = grid.variables['clon']
clat = grid.variables['clat']
ncells_grid = len(clon)

neighbors = grid.variables['neighbor_cell_index']

qc_fill = np.zeros((ncells_grid))

neighbors = np.ma.masked_where(neighbors==-1,neighbors)

print min(clon)
print min(clat)
print max(clon)
print max(clat)


for i in range(0,ncells_grid):
#for i in range(0,10):
    #print 'i:',i
    idx_neighbors = neighbors[:,i]
    #print idx_neighbors
    #print qc_fill[idx_neighbors[0]]
    #print qc_fill[idx_neighbors[1]]
    #print qc_fill[idx_neighbors[2]]
    qc_sum = sum(qc_fill[idx_neighbors])
    if qc_sum<1.:
        qc_fill[i+1] = 1.
     #   print 'in if loop:',i

#qc_fill[0:11] = 1

qc_real = np.zeros((ncells_grid))
for j in range(0,ncells_grid):
    if qc_fill[j]==1:
        qc_real[j-1] = 1




subsub = NetCDFFile('subsub_regular.nc','w')

ncells = subsub.createDimension('ncells',ncells_grid)
time = subsub.createDimension('time',0)
height = subsub.createDimension('height',0)

qc = subsub.createVariable('qc',np.float32,['time','height','ncells'])

qc[0,0,:] = qc_real
qc.units = 'kg kg-1'
qc.long_name = 'specific_cloud_water_content'

subsub.close()
