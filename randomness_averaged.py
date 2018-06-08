#!/usr/bin/python

# Find minimal distances between clouds in one bin, average these per bin
# Compute geometric and arithmetical mean between all clouds per bin

from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
import numpy as np
from numpy import ma
from distance_methods import distances,SCAI,COP, NNCDF
from Randomfield import randomfield

cusize = NetCDFFile(
         '/home/vanlaar/HDCP2data/TA_dom4/python_output/cusize_output.nc')

size = cusize.variables['size']

grid = NetCDFFile('../NarvalDom2_NestE-R02B14_DOM03.nc')
clon_grid = grid.variables['clon']
clat_grid = grid.variables['clat']

Nmax = len(clon_grid)/2
lengthscale = 440000

ntime = len(cusize.dimensions['time'])
print ntime


#######################################################


D0_all = np.zeros((ntime,len(size)))
D1_all = np.zeros((ntime,len(size)))
nclouds_bin_all = np.zeros([ntime,len(size)])
nncdf_real = np.zeros([ntime,99])
bins_nncdf = np.zeros([ntime,100])
nclouds_all = np.zeros(ntime)
SCAI_0_all = np.zeros(ntime)
SCAI_1_all = np.zeros(ntime)
COP_all = np.zeros(ntime)


#for time in range(0,ntime):
for time in range(40,41):
    print 'time:',time

    ncloudsint = cusize.variables['nclouds'][time]
    nclouds_real = np.real(ncloudsint)

    print ncloudsint

    cloudlon = cusize.variables['cloud_lon'][time,0:ncloudsint]
    cloudlat = cusize.variables['cloud_lat'][time,0:ncloudsint]
    cloudsize = cusize.variables['cloud_size'][time,0:ncloudsint]
    cloudbin = np.ceil(cloudsize/size[0])
    CSD = cusize.variables['CSD'][0,:]
    CSD_normalized = CSD/nclouds_real

    filledbin = int(max(cloudbin))
    
    output_distances = distances(filledbin,cloudlon,cloudlat,cloudbin,size,ncloudsint)

    D0_all[time] = output_distances[0]
    D1_all[time] = output_distances[1]
    nclouds_bin_all[time] = output_distances[4]
    nclouds_all[time] = nclouds_real

    SCAIs = SCAI(cloudlon,cloudlat,ncloudsint,Nmax,lengthscale)
    COP_output = COP(cloudlon,cloudlat,cloudsize)

    SCAI_0_all[time] = SCAIs[0]
    SCAI_1_all[time] = SCAIs[1]
    COP_all[time] = COP_output

    print 'Real field:'
    print 'SCAI 0:', SCAI_0_all[time]
    print 'SCAI 1:', SCAI_1_all[time]
    print 'COP:', COP_all[time]
    
    nncdf_output = NNCDF(cloudlon,cloudlat)
    nncdf_real[time,:] = nncdf_output[0]
    bins_nncdf[time,:] = nncdf_output[1]
   

    #####Random#####

    #nrbins = len(size)
    #probability = nclouds_bin_all[time-41]/nclouds_real # probability of drawing from bin

    #randomfield_output = randomfield(ncloudsint,nrbins,clon_grid,clat_grid,probability)
    #cloud_lon_random = randomfield_output[0]
    #cloud_lat_random = randomfield_output[1]
    #cloud_bin_random = np.array(randomfield_output[2])
    
    #distances_random = distances(filledbin,cloud_lon_random,cloud_lat_random,cloud_bin_random,size,nclouds_av)
   
    #lon_random_min = cloud_lon_random_min[time-41,0:ncloudsint]
    #lat_random_min = cloud_lat_random_min[time-41,0:ncloudsint]
    #bin_random_min = cloud_bin_random_min[time-41,0:ncloudsint]

    #lon_random_max = cloud_lon_random_max[time-41,0:ncloudsint]
    #lat_random_max = cloud_lat_random_max[time-41,0:ncloudsint]
    #bin_random_max = cloud_bin_random_max[time-41,0:ncloudsint]
 
    #SCAIs_min = SCAI(lon_random_min,lat_random_min,ncloudsint,Nmax,lengthscale)
    #SCAI_0_random_all_min[time-41] = SCAIs_min[0]
    #SCAI_1_random_all_min[time-41] = SCAIs_min[1]

    #SCAIs_max = SCAI(lon_random_min,lat_random_min,ncloudsint,Nmax,lengthscale)
    #SCAI_0_random_all_max[time-41] = SCAIs_max[0]
    #SCAI_1_random_all_max[time-41] = SCAIs_max[1]
    
    #cloud_size_random = bin_random_min*size[0]
    #COP_random_output_min = COP(lon_random_min,lat_random_min,cloud_size_random)
    #COP_random_all_min[time-41] = COP_random_output_min

    #cloud_size_random = bin_random_max*size[0]
    #COP_random_output_max = COP(lon_random_max,lat_random_max,cloud_size_random)
    
    #COP_random_all_max[time-41] = COP_random_output_max
    
    #print 'Random:'
    #print 'SCAI 0:', SCAI_0_random_all_min[time-41]
    #print 'SCAI 1:', SCAI_1_random_all_min[time-41]
    #print 'COP:', COP_random_all_min[time-41]

    #nncdf_random_min = NNCDF(lon_random_min, lat_random_min, size)
    #area_min[time-41] = np.trapz(nncdf_real, nncdf_random_min)
    #print area_min[time-41]
    #nncdf_random_max = NNCDF(lon_random_max, lat_random_max, size)
    #area_max[time-41] = np.trapz(nncdf_real, nncdf_random_max)
    #print area_max[time-41]


#D0 = np.mean(D0_all,axis=0)
#D1 = np.mean(D1_all,axis=0)
#nclouds_av = int(np.mean(nclouds_all))
#nclouds_bin_av = np.mean(nclouds_bin_all,axis=0)


#filledbin = np.argmin(mindistance_mean)
#slope, intercept, r_value, p_value, std_err = stats.linregress(size[0:filledbin],mindistance_mean[0:filledbin])
#line = intercept + slope*size

#############################################################
#### Randomness ####
############################################################

#unique, counts = np.unique(cloud_bin_random, return_counts=True)
#frequency_random = np.array(counts).astype(np.float64)
#frequency_randomT = frequency_random[:]/float(nclouds_real)
#
#hn_random = np.zeros(nrbins,float)
#hn_random[unique-1] = frequency_randomT[:]
#print sum(hn_random)
#hn_random_log = hn_random
#hn_random_log[hn_random>0] = np.log10(hn_random[hn_random>0])
#
#cloud_bin_real = nclouds_bin_av.astype(np.float64)
#frequencyT = cloud_bin_real/float(nclouds_real)
#print sum(frequencyT)
#hn_real_log = frequencyT
#hn_real_log[frequencyT>0] = np.log10(frequencyT[frequencyT>0])
#
#size_log = np.log10(size)
#
#plt.figure(figsize=(10,8))
#plt.axis([2, 4.5, -4, -0.1])
#plt.scatter(size_log,hn_real_log,c='k')
#plt.scatter(size_log,hn_random_log,c='b')
#plt.savefig('Figures/CSD_realrandom.pdf')




###############################################################
###### Netcdffile#######
######################################################

random_parameters = NetCDFFile('random_averaged_test.nc','w')

time = random_parameters.createDimension('time',ntime)
size = random_parameters.createDimension('size',len(size))
nrbins_nncdf = random_parameters.createDimension('nrbins_nncdf',99)
binedges = random_parameters.createDimension('binedges',100)

SCAI_0_real = random_parameters.createVariable('SCAI_0_real',np.float64,['time'])
SCAI_0_real[:] = SCAI_0_all[:]

SCAI_1_real = random_parameters.createVariable('SCAI_1_real',np.float64,['time'])
SCAI_1_real[:] = SCAI_1_all[:] 

COP_real = random_parameters.createVariable('COP_real',np.float64,['time'])
COP_real[:] = COP_all[:]

nclouds = random_parameters.createVariable('nclouds',np.float64,['time'])
nclouds[:] = nclouds_all[:]

nncdf = random_parameters.createVariable('nncdf',np.float64,['time','nrbins_nncdf'])
nncdf[:] = nncdf_real[:]

nncdf_binedges = random_parameters.createVariable('nncdf_binedges',np.float64,['time','binedges'])
nncdf_binedges[:] = bins_nncdf[:]


random_parameters.close()

