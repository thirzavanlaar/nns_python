#!/usr/bin/python


from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.tri as tri
from scipy import spatial

cusize = NetCDFFile('/home/vanlaar/HDCP2data/cusize_subsub_regular.nc')
grid = NetCDFFile('subsubgrid.nc')
inputfile = NetCDFFile('subsub_regular.nc')

clon = np.rad2deg(grid.variables['clon'])
clat = np.rad2deg(grid.variables['clat'])

#vlon = np.rad2deg(grid.variables['clon_bnds'])
#vlat = np.rad2deg(grid.variables['clat_bnds'])

x = clon
y = clat

x_y = np.vstack((clon,clat)).T
Delau = spatial.Delaunay(x_y,furthest_site=True)

print Delau.simplices.shape

#triang = tri.Triangulation(x,y,triangles=tri.simplices)
triang = tri.Triangulation(x,y,triangles=Delau.simplices)


#plt.axis([-56.975, -56.95, 13.36, 13.38])
#plt.triplot(x_y[:,0],x_y[:,1],tri.simplices.copy())
#plt.plot(x_y[:,0], x_y[:,1], 'o')
#plt.show()



qc_all = inputfile.variables['qc']
qc = qc_all[0,0,:]

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

plt.figure(figsize=(10,8))
#plt.axis([-60, -55.5, 12.3, 14])
#plt.xlabel('Cloud size [m]')
#plt.ylabel('Nearest-neighbour distance [m]')
#plt.fill_between(size,mindistance_plus,mindistance_minus,alpha=0.3,color='red')
#plt.scatter(cloud_lon_deg,cloud_lat_deg,s=cloud_bin)
#plt.tricontourf(triang,qc,2,c='blue')
plt.tricontourf(x,y,triang,qc)
#plt.plot(size,line,color='black')
plt.savefig('../Figures/projected_clouds_subsub_regular.png')










