#!/usr/bin/python


from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
import numpy as np
from distance_methods import distances,NNCDF


random_overlap = NetCDFFile('random_overlap.nc')

cloud_bin_overlap = random_overlap.variables['cloud_bin_overlap'][:]
cloud_lon_overlap = random_overlap.variables['cloud_lon_overlap'][:]
cloud_lat_overlap = random_overlap.variables['cloud_lat_overlap'][:]
overlap_var = random_overlap.variables['overlap_var'][:]

cloud_lon_real = random_overlap.variables['cloud_lon_real'][:]
cloud_lat_real = random_overlap.variables['cloud_lat_real'][:]
size = random_overlap.variables['size'][:]

max_idx = np.argmax(overlap_var)
min_idx = np.argmin(overlap_var)

print overlap_var[max_idx]
print overlap_var[min_idx]

nncdf_real = NNCDF(cloud_lon_real,cloud_lat_real,size)

nncdf_random_min = NNCDF(cloud_lon_overlap[min_idx,:],cloud_lat_overlap[min_idx,:],size)
nncdf_random_max = NNCDF(cloud_lon_overlap[max_idx,:],cloud_lat_overlap[max_idx,:],size)





####################################################################
#### Figures ####
####################################################################

orange = (1.,0.38,0.01)
blue = (0.53,0.81,1)


plt.figure(figsize=(10,8))
#plt.axis([50000,220000, 50000, 220000])
plt.xlabel('nncdf random', fontsize=15)
plt.ylabel('nncdf real', fontsize=15)
plt.scatter(nncdf_random_max,nncdf_real, label='max')
plt.scatter(nncdf_random_min,nncdf_real, label='min')
plt.plot(nncdf_random_min,nncdf_random_min,c='black')
plt.legend()
plt.savefig('Figures/nncdf_maxmin.pdf')


