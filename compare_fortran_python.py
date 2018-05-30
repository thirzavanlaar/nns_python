#!/usr/bin/python

# apply hierarchical clustering to output cusize based on the 'Interaction Potential'

from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, centroid
import numpy as np
from scipy.spatial import distance
from haversine import haversine

cusize = NetCDFFile('/home/vanlaar/HDCP2data/TA_dom4/python_output/cusize_output_time44_square.nc')

# Import data from netcdf file:
cloudlon_p = cusize.variables['cloud_lon'][:]
cloudlat_p = cusize.variables['cloud_lat'][:]
nclouds_cusize_p  = cusize.variables['nclouds']
cloud_size_p = cusize.variables['cloud_size'][:]
size_p = cusize.variables['size'][:]
nclouds_p = int(nclouds_cusize_p[0])
CSD = cusize.variables['CSD'][:]


# Adjust data format for later use:
#cloud_lon_p = cloudlon_p[0,0:nclouds_p]
#cloud_lat_p = cloudlat_p[0,0:nclouds_p]
#cloud_size_p = cloud_bin_p[0,:nclouds_p]*size_p[0]



cusize = NetCDFFile('/home/vanlaar/HDCP2data/TA_dom4/cusize_output_time44.nc')

# Import data from netcdf file:
cloudlon_f = cusize.variables['cloud_lon']
cloudlat_f = cusize.variables['cloud_lat']
nclouds_cusize_f  = cusize.variables['nclouds']
cloud_bin_f = cusize.variables['cloud_bin']
size_f = cusize.variables['size']
nclouds_f = int(nclouds_cusize_f[0])
hn = cusize.variables['hn'][:]

# Adjust data format for later use:
cloud_lon_f = cloudlon_f[0,0:nclouds_f]
cloud_lat_f = cloudlat_f[0,0:nclouds_f]
cloud_size_f = cloud_bin_f[0,:nclouds_f]*size_f[0]


print 'nclouds python:', nclouds_p
print 'nclouds fortran:', nclouds_f

CSD_log = np.where(CSD>0, np.log10(CSD),0.)
hn_log = np.where(hn>0,np.log10(hn/nclouds_cusize_f),0.)
size_f_log = np.where(size_f>0, np.log10(size_f),0.)
size_p_log = np.where(size_p>0, np.log10(size_p),0.)



plt.figure()
axes = plt.gca()
axes.set_xlim([1.5,4.5])
axes.set_ylim([-6.5,-2])
#axes.set_yscale('log')
#axes.set_xscale('log')
plt.scatter(size_f_log,hn_log, label = 'Fortran')
plt.scatter(size_p_log,CSD_log, label = 'Python')
plt.legend()
plt.title('Clouds are a square')
plt.ylabel('log(N*(l))')
plt.xlabel('log(l)')
plt.savefig('Figures/compare_p_f_square.pdf')
