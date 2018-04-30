#!/usr/bin/python


from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
import numpy as np
from numpy import ma
from scipy import stats 
from haversine import haversine
from scipy.spatial import distance
from distance_methods import distances
import matplotlib.tri as tri

cusize = NetCDFFile(
         '/home/vanlaar/HDCP2data/TA_dom4/cusize_output_time41.nc')

size = cusize.variables['size']

begin_time = 41
end_time = 48


cloudlon = cusize.variables['cloud_lon']
cloudlat = cusize.variables['cloud_lat']
size = cusize.variables['size']
cloud_bin = cusize.variables['cloud_bin']
nclouds_cusize = cusize.variables['nclouds']
hn = cusize.variables['hn']

nclouds = int(nclouds_cusize[0])
cloud_lon = cloudlon[0,0:nclouds]
cloud_lat = cloudlat[0,0:nclouds]
filledbin=np.argmin(hn[0,:])  # last bin with clouds, rest is empty
cloud_size = cloud_bin*size[0]

cloud_lon_deg = np.rad2deg(cloud_lon)
cloud_lat_deg = np.rad2deg(cloud_lat)

print 'cusize:'
print min(cloud_lon)
print max(cloud_lon)

plt.figure(figsize=(10,8))
#plt.axis([-60, -55.5, 12.3, 14])
#plt.xlabel('Cloud size [m]')
#plt.ylabel('Nearest-neighbour distance [m]')
#plt.fill_between(size,mindistance_plus,mindistance_minus,alpha=0.3,color='red')
plt.scatter(cloud_lon_deg,cloud_lat_deg,s=cloud_bin)
#plt.plot(size,line,color='black')
plt.savefig('Figures/projected_clouds_TA_cusize.png')


# now directly from ICON output and not from cusize

ICON = NetCDFFile('/data/inscape/vanlaar/HDCP2_data/input/TA_dom4_qc_time41.nc')

clon_rad = ICON.variables['clon']

clon = np.rad2deg(ICON.variables['clon'])
clat = np.rad2deg(ICON.variables['clat'])

vlon = np.rad2deg(ICON.variables['clon_bnds'])
vlat = np.rad2deg(ICON.variables['clat_bnds'])

x = clon
y = clat

triang = tri.Triangulation(x,y)

qc_output = ICON.variables['qc']  # [time,lev,cell]
qc = qc_output[0,:,:]

column_sum = np.sum(qc,axis=0)
qc_mask = np.where(column_sum>0.0,1,0)

#print 'ICON:'
#print min(clon_rad)
#print max(clon_rad)

plt.figure(figsize=(10,8))
#plt.axis([-60, -55.5, 12.3, 14])
#plt.xlabel('Cloud size [m]')
#plt.ylabel('Nearest-neighbour distance [m]')
#plt.fill_between(size,mindistance_plus,mindistance_minus,alpha=0.3,color='red')
plt.tricontourf(triang,qc_mask,c='blue')
#plt.plot(size,line,color='black')
plt.savefig('Figures/projected_clouds_TA_ICON.png')










