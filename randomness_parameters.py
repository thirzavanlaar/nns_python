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
         '/home/vanlaar/HDCP2data/TA_dom4/cusize_output_time41.nc')

size = cusize.variables['size']

grid = NetCDFFile('../NarvalDom2_NestE-R02B14_DOM03.nc')
clon_grid = grid.variables['clon']
clat_grid = grid.variables['clat']

Nmax = len(clon_grid)/2
lengthscale = 440000

begin_time = 41
end_time = 48

# Read random min and max overlap fields:

random = NetCDFFile('random_overlap.nc')

cloud_bin_random_min = random.variables['cloud_bin_overlap_min'][:]
cloud_lon_random_min = random.variables['cloud_lon_overlap_min'][:]
cloud_lat_random_min = random.variables['cloud_lat_overlap_min'][:]
cloud_bin_random_max = random.variables['cloud_bin_overlap_max'][:]
cloud_lon_random_max = random.variables['cloud_lon_overlap_max'][:]
cloud_lat_random_max = random.variables['cloud_lat_overlap_max'][:]

print cloud_lon_random_min.shape


#######################################################


D0_all = np.zeros((end_time-begin_time+1,len(size)))
D1_all = np.zeros((end_time-begin_time+1,len(size)))
nclouds_bin_all = np.zeros((end_time-begin_time+1,len(size)))
nclouds_all = np.zeros((end_time-begin_time+1))
SCAI_0_all = np.zeros((end_time-begin_time+1))
SCAI_1_all = np.zeros((end_time-begin_time+1))
COP_all = np.zeros((end_time-begin_time+1))
COP_random_all_min = np.zeros((end_time-begin_time+1))
SCAI_0_random_all_min = np.zeros((end_time-begin_time+1))
SCAI_1_random_all_min = np.zeros((end_time-begin_time+1))
area_min = np.zeros((end_time-begin_time+1))
COP_random_all_max = np.zeros((end_time-begin_time+1))
SCAI_0_random_all_max = np.zeros((end_time-begin_time+1))
SCAI_1_random_all_max = np.zeros((end_time-begin_time+1))
area_max = np.zeros((end_time-begin_time+1))


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
    cloud_lon = cloudlon[0,0:ncloudsint]
    cloud_lat = cloudlat[0,0:ncloudsint]
    cloud_bin = cloudbin[0,0:ncloudsint] 
    filledbin = int(max(cloud_bin))
    
    output_distances = distances(filledbin,cloud_lon,cloud_lat,cloud_bin,size,ncloudsint)

    D0_all[time-41] = output_distances[0]
    D1_all[time-41] = output_distances[1]
    nclouds_bin_all[time-41] = output_distances[4]
    nclouds_all[time-41] = nclouds_real

    SCAIs = SCAI(cloud_lon,cloud_lat,ncloudsint,Nmax,lengthscale)

    cloud_size = cloud_bin*size[0]
    COP_output = COP(cloud_lon,cloud_lat,cloud_size)

    SCAI_0_all[time-41] = SCAIs[0]
    SCAI_1_all[time-41] = SCAIs[1]
    COP_all[time-41] = COP_output

    print 'Real field:'
    print 'SCAI 0:', SCAI_0_all[time-41]
    print 'SCAI 1:', SCAI_1_all[time-41]
    print 'COP:', COP_all[time-41]

    nncdf_real = NNCDF(cloud_lon,cloud_lat,size)

    #####Random#####

    nrbins = len(size)
    #probability = nclouds_bin_all[time-41]/nclouds_real # probability of drawing from bin

    #randomfield_output = randomfield(ncloudsint,nrbins,clon_grid,clat_grid,probability)
    #cloud_lon_random = randomfield_output[0]
    #cloud_lat_random = randomfield_output[1]
    #cloud_bin_random = np.array(randomfield_output[2])
    
    #distances_random = distances(filledbin,cloud_lon_random,cloud_lat_random,cloud_bin_random,size,nclouds_av)
   
    lon_random_min = cloud_lon_random_min[time-41,0:ncloudsint]
    lat_random_min = cloud_lat_random_min[time-41,0:ncloudsint]
    bin_random_min = cloud_bin_random_min[time-41,0:ncloudsint]

    lon_random_max = cloud_lon_random_max[time-41,0:ncloudsint]
    lat_random_max = cloud_lat_random_max[time-41,0:ncloudsint]
    bin_random_max = cloud_bin_random_max[time-41,0:ncloudsint]
 
    SCAIs_min = SCAI(lon_random_min,lat_random_min,ncloudsint,Nmax,lengthscale)
    SCAI_0_random_all_min[time-41] = SCAIs_min[0]
    SCAI_1_random_all_min[time-41] = SCAIs_min[1]

    SCAIs_max = SCAI(lon_random_min,lat_random_min,ncloudsint,Nmax,lengthscale)
    SCAI_0_random_all_max[time-41] = SCAIs_max[0]
    SCAI_1_random_all_max[time-41] = SCAIs_max[1]
    
    cloud_size_random = bin_random_min*size[0]
    COP_random_output_min = COP(lon_random_min,lat_random_min,cloud_size_random)
    COP_random_all_min[time-41] = COP_random_output_min

    cloud_size_random = bin_random_max*size[0]
    COP_random_output_max = COP(lon_random_max,lat_random_max,cloud_size_random)
    
    COP_random_all_max[time-41] = COP_random_output_max
    
    print 'Random:'
    print 'SCAI 0:', SCAI_0_random_all_min[time-41]
    print 'SCAI 1:', SCAI_1_random_all_min[time-41]
    print 'COP:', COP_random_all_min[time-41]

    nncdf_random_min = NNCDF(lon_random_min, lat_random_min, size)
    area_min[time-41] = np.trapz(nncdf_real, nncdf_random_min)
    print area_min[time-41]
    nncdf_random_max = NNCDF(lon_random_max, lat_random_max, size)
    area_max[time-41] = np.trapz(nncdf_real, nncdf_random_max)
    print area_max[time-41]


D0 = np.mean(D0_all,axis=0)
D1 = np.mean(D1_all,axis=0)
nclouds_av = int(np.mean(nclouds_all))
nclouds_bin_av = np.mean(nclouds_bin_all,axis=0)


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


####################################################################
#### Figures ####
####################################################################

orange = (1.,0.38,0.01)
blue = (0.53,0.81,1)


plt.figure(figsize=(10,8))
#plt.axis([50000,220000, 50000, 220000])
plt.xlabel('D1')
plt.ylabel('D0')
plt.scatter(D1,D0,color='k')
plt.savefig('Figures/D1-D0.pdf')


#plt.figure(figsize=(10,8))
#plt.scatter(cloud_lon_random,cloud_lat_random,s=cloud_bin_random)
#plt.savefig('randomfield.pdf')


###############################################################
###### Netcdffile#######
######################################################

random_parameters = NetCDFFile('random_parameters.nc','w')

time = random_parameters.createDimension('time',end_time-begin_time+1)


SCAI_0_real = random_parameters.createVariable('SCAI_0_real',np.float64,['time'])
SCAI_0_real[:] = SCAI_0_all[:]

SCAI_1_real = random_parameters.createVariable('SCAI_1_real',np.float64,['time'])
SCAI_1_real[:] = SCAI_1_all[:] 

COP_real = random_parameters.createVariable('COP_real',np.float64,['time'])
COP_real[:] = COP_all[:]

SCAI_0_random_min = random_parameters.createVariable('SCAI_0_random_min',np.float64,['time'])
SCAI_0_random_min[:] = SCAI_0_random_all_min[:]

SCAI_1_random_min = random_parameters.createVariable('SCAI_1_random_min',np.float64,['time'])
SCAI_1_random_min[:] = SCAI_1_random_all_min[:]

COP_random_min = random_parameters.createVariable('COP_random_min',np.float64,['time'])
COP_random_min[:] = COP_random_all_min[:]

SCAI_0_random_max = random_parameters.createVariable('SCAI_0_random_max',np.float64,['time'])
SCAI_0_random_max[:] = SCAI_0_random_all_max[:]

SCAI_1_random_max = random_parameters.createVariable('SCAI_1_random_max',np.float64,['time'])
SCAI_1_random_max[:] = SCAI_1_random_all_max[:]

COP_random_max = random_parameters.createVariable('COP_random_max',np.float64,['time'])
COP_random_max[:] = COP_random_all_max[:]

nclouds = random_parameters.createVariable('nclouds',np.float64,['time'])
nclouds[:] = nclouds_all[:]

Iorg_min = random_parameters.createVariable('Iorg_min',np.float64,['time'])
Iorg_min[:] = area_min[:]
Iorg_max = random_parameters.createVariable('Iorg_max',np.float64,['time'])
Iorg_max[:] = area_max[:]

random_parameters.close()

