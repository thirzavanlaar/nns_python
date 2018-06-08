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
         '/home/vanlaar/HDCP2data/TA_dom4/python_output/cusize_output.nc')

size = cusize.variables['size']

grid = NetCDFFile('../NarvalDom2_NestE-R02B14_DOM03.nc')
clon_grid = grid.variables['clon']
clat_grid = grid.variables['clat']

Nmax = len(clon_grid)/2
lengthscale = 440000

ntime = len(cusize.dimensions['time'])


binwidth = size[0]
area = np.zeros((ntime))
cloud_lon_min = np.zeros((ntime,Nmax))
cloud_lat_min = np.zeros((ntime,Nmax))
cloud_bin_min = np.zeros((ntime,Nmax))
random_overlap_output = np.zeros((ntime))
cloud_lon_max = np.zeros((ntime,Nmax))
cloud_lat_max = np.zeros((ntime,Nmax))
cloud_bin_max = np.zeros((ntime,Nmax))


for time in range(0,ntime):
#for time in range(0,1):
    print 'time:',time

        
    ncloudsint = cusize.variables['nclouds'][time]
    nclouds_real = np.real(ncloudsint)

    cloudlon = cusize.variables['cloud_lon'][time,0:ncloudsint]
    cloudlat = cusize.variables['cloud_lat'][time,0:ncloudsint]
    cloudsize = cusize.variables['cloud_size'][time,0:ncloudsint]
    cloudbin = np.ceil(cloudsize/size[0])

    filledbin = int(max(cloudbin))
    
    output_distances = distances(filledbin,cloudlon,cloudlat,cloudbin,size,ncloudsint)

    nclouds_bin_all = output_distances[4]
    nclouds_all = nclouds_real
   
    print nclouds_bin_all.shape
    print nclouds_all.shape
 
    ####################
    ####Random field####
    ####################

    #probability = nclouds_bin_all[time-41]/nclouds_real
    probability = nclouds_bin_all/nclouds_real
    nrbins = len(size)

    print probability.shape

    overlap_output = overlap(ncloudsint,nrbins,clon_grid,clat_grid,probability,binwidth)
    #randomfield_output = randomfield(ncloudsint,nrbins,clon_grid,clat_grid,probability)
   
    cloud_lon_random = overlap_output[0]
    cloud_lat_random = overlap_output[1]
    cloud_bin_random = np.array(overlap_output[2])
    random_overlap_output = np.array(overlap_output[3])

    max_idx = np.argmax(random_overlap_output)
    min_idx = np.argmin(random_overlap_output)

    cloud_lon_min[time,0:ncloudsint] = cloud_lon_random[min_idx,:]
    cloud_lat_min[time,0:ncloudsint] = cloud_lat_random[min_idx,:]
    cloud_bin_min[time,0:ncloudsint] = cloud_bin_random[min_idx,:]

    cloud_lon_max[time,0:ncloudsint] = cloud_lon_random[max_idx,:]
    cloud_lat_max[time,0:ncloudsint] = cloud_lat_random[max_idx,:]
    cloud_bin_max[time,0:ncloudsint] = cloud_bin_random[max_idx,:]


    
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
nr_iterations = random_overlap.createDimension('nr_iterations',ntime)
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


