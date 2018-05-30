#!/usr/bin/python

# Find minimal distances between clouds in one bin, average these per bin
# Compute geometric and arithmetical mean between all clouds per bin

from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
import numpy as np
from numpy import ma
from distance_methods import distances,NNCDF
from Randomfield import randomfield, randomfield_nooverlap, overlap

cusize = NetCDFFile(
         '/home/vanlaar/HDCP2data/TA_dom4/cusize_output_time41.nc')

size = cusize.variables['size']

grid = NetCDFFile('../NarvalDom2_NestE-R02B14_DOM03.nc')
clon_grid = grid.variables['clon']
clat_grid = grid.variables['clat']

Nmax = len(clon_grid)/2
lengthscale = 440000

begin_time = 41
end_time = 48

binwidth = size[0]
D0_all = np.zeros((end_time-begin_time+1,len(size)))
D1_all = np.zeros((end_time-begin_time+1,len(size)))
nclouds_bin_all = np.zeros((end_time-begin_time+1,len(size)))
nclouds_all = np.zeros((end_time-begin_time+1))
SCAI_0_all = np.zeros((end_time-begin_time+1))
SCAI_1_all = np.zeros((end_time-begin_time+1))
area = np.zeros((end_time-begin_time+1))
cloud_lon_min = np.zeros((end_time-begin_time+1,Nmax))
cloud_lat_min = np.zeros((end_time-begin_time+1,Nmax))
cloud_bin_min = np.zeros((end_time-begin_time+1,Nmax))
random_overlap_output = np.zeros((end_time-begin_time+1))
cloud_lon_max = np.zeros((end_time-begin_time+1,Nmax))
cloud_lat_max = np.zeros((end_time-begin_time+1,Nmax))
cloud_bin_max = np.zeros((end_time-begin_time+1,Nmax))
cloud_lon = np.zeros((end_time-begin_time+1,Nmax))
cloud_lat = np.zeros((end_time-begin_time+1,Nmax))
cloud_bin = np.zeros((end_time-begin_time+1,Nmax))

for time in range(begin_time,end_time+1):
    print 'time:',time

    cusize = NetCDFFile(
         '/home/vanlaar/HDCP2data/TA_dom4/cusize_output_time'+str(time)+'.nc')

    cloudlon = cusize.variables['cloud_lon']
    cloudlat = cusize.variables['cloud_lat']
    nclouds_real = cusize.variables['nclouds'][0]
    size = cusize.variables['size']
    cloudbin = cusize.variables['cloud_bin']
    hn = cusize.variables['hn']
    hn_normalized_loop = hn/nclouds_real
    ncloud_bin = cusize.variables['ncloud_bin'][:]

    ncloudsint = int(nclouds_real)
    cloud_lon[time-41,0:ncloudsint] = cloudlon[0,0:ncloudsint]
    cloud_lat[time-41,0:ncloudsint] = cloudlat[0,0:ncloudsint]
    cloud_bin[time-41,0:ncloudsint] = cloudbin[0,0:ncloudsint] 
    filledbin = int(max(cloud_bin[time-41,0:ncloudsint]))
    
    output_distances = distances(filledbin,cloud_lon[time-41,0:ncloudsint],cloud_lat[time-41,0:ncloudsint],cloud_bin[time-41,0:ncloudsint],size,ncloudsint)

    #nclouds_bin_all[time-41] = output_distances[4]
    #nclouds_all[time-41] = nclouds_real
    nclouds_bin_all[0] = output_distances[4]
    nclouds_all[0] = nclouds_real
    
    #nncdf_real = NNCDF(cloud_lon,cloud_lat,size)

    ####################
    ####Random field####
    ####################

    #probability = nclouds_bin_all[time-41]/nclouds_real
    probability = nclouds_bin_all[0]/nclouds_real
    nrbins = len(size)

    binwidth = size[0]

    overlap_output = overlap(ncloudsint,nrbins,clon_grid,clat_grid,probability,binwidth)
    #randomfield_output = randomfield(ncloudsint,nrbins,clon_grid,clat_grid,probability)
   
    cloud_lon_random = overlap_output[0]
    cloud_lat_random = overlap_output[1]
    cloud_bin_random = np.array(overlap_output[2])
    random_overlap_output = np.array(overlap_output[3])

    max_idx = np.argmax(random_overlap_output)
    min_idx = np.argmin(random_overlap_output)

    cloud_lon_min[time-41,0:ncloudsint] = cloud_lon_random[min_idx,:]
    cloud_lat_min[time-41,0:ncloudsint] = cloud_lat_random[min_idx,:]
    cloud_bin_min[time-41,0:ncloudsint] = cloud_bin_random[min_idx,:]

    cloud_lon_max[time-41,0:ncloudsint] = cloud_lon_random[max_idx,:]
    cloud_lat_max[time-41,0:ncloudsint] = cloud_lat_random[max_idx,:]
    cloud_bin_max[time-41,0:ncloudsint] = cloud_bin_random[max_idx,:]


    
    #distances_random = distances(filledbin,cloud_lon_random,cloud_lat_random,cloud_bin_random,size,nrclouds)
    
    #nncdf_random = NNCDF(cloud_lon_random,cloud_lat_random,size)

    #area[time-41] = np.trapz(nncdf_real, nncdf_random)
    #area[0] = np.trapz(nncdf_real, nncdf_random)
    #print area
    



#D0 = np.mean(D0_all,axis=0)
#D1 = np.mean(D1_all,axis=0)
#nclouds_av = int(np.mean(nclouds_all))
#nclouds_bin_av = np.mean(nclouds_bin_all,axis=0)




###########################################################
###Netcdffile#############
###############################################

random_overlap = NetCDFFile('random_overlap.nc','w')

#nr_iterations = random_overlap.createDimension('nr_iterations',100)
nr_iterations = random_overlap.createDimension('nr_iterations',end_time-begin_time+1)
nclouds = random_overlap.createDimension('nclouds',Nmax)
size_dim = random_overlap.createDimension('size',0)

cloud_bin_overlap_min = random_overlap.createVariable('cloud_bin_overlap_min',np.float64,['nr_iterations','nclouds'])
cloud_bin_overlap_min[:] = cloud_bin_min[:]

cloud_lon_overlap_min = random_overlap.createVariable('cloud_lon_overlap_min',np.float64,['nr_iterations','nclouds'])
cloud_lon_overlap_min[:] = cloud_lon_min[:]

cloud_lat_overlap_min = random_overlap.createVariable('cloud_lat_overlap_min',np.float64,['nr_iterations','nclouds'])
cloud_lat_overlap_min[:] = cloud_lat_min[:]

cloud_bin_overlap_max = random_overlap.createVariable('cloud_bin_overlap_max',np.float64,['nr_iterations','nclouds'])
cloud_bin_overlap_max[:] = cloud_bin_max[:]

cloud_lon_overlap_max = random_overlap.createVariable('cloud_lon_overlap_max',np.float64,['nr_iterations','nclouds'])
cloud_lon_overlap_max[:] = cloud_lon_max[:]

cloud_lat_overlap_max = random_overlap.createVariable('cloud_lat_overlap_max',np.float64,['nr_iterations','nclouds'])
cloud_lat_overlap_max[:] = cloud_lat_max[:]

cloud_lon_real = random_overlap.createVariable('cloud_lon_real',np.float64,['nr_iterations','nclouds'])
cloud_lon_real[:] = cloud_lon[:]

cloud_lat_real = random_overlap.createVariable('cloud_lat_real',np.float64,['nr_iterations','nclouds'])
cloud_lat_real[:] = cloud_lat[:]

cloud_bin_real = random_overlap.createVariable('cloud_bin_real',np.float64,['nr_iterations','nclouds'])
cloud_bin_real[:] = cloud_bin[:]

size_var = random_overlap.createVariable('size',np.float64,['size'])
size_var[:] = size[:]



random_overlap.close()

####################################################################
#### Figures ####
####################################################################

orange = (1.,0.38,0.01)
blue = (0.53,0.81,1)


#plt.figure(figsize=(10,8))
#plt.axis([50000,220000, 50000, 220000])
#plt.xlabel('D1')
#plt.ylabel('D0')
#plt.scatter(D1,D0,color='k')
#plt.savefig('Figures/D1-D0.pdf')


