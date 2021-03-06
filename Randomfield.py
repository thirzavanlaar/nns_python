#!/usr/bin/python

from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
import numpy as np
import random
from haversine import haversine
from scipy.spatial import distance
import itertools

def randomfield(nclouds,nrbins,clon,clat,probability):

    """
    Compute a random cloud field

    Clouds are sampled from the CSD of the real field (represented by 'probability' and placed randomly in the field (clon and clat). Outputs are the coordinates of the clouds and their size (represented by the binnr they are in).
    """

    ncells_grid = len(clon)
    ncells_range = np.arange(0,ncells_grid,1)
    random_cells = random.sample(ncells_range, nclouds)

    cloud_lon_random = clon[random_cells]
    cloud_lat_random = clat[random_cells]

    nbins_range = np.arange(0,nrbins,1)
    cloud_bin_random = np.random.choice(nbins_range,nclouds, True, probability)
    cloud_bin_random = cloud_bin_random +1

    return cloud_lon_random,cloud_lat_random,cloud_bin_random



def randomfield_nooverlap(nclouds,nrbins,clon,clat,probability,binwidth):

    """
    Compute a random cloud field without overlap

    Clouds are sampled from the CSD of the real field and randomly placed in the domain. See also 'randomfield'. When a field is generated it is checked for overlap. When there is overlap, a new one is made. This continues as long as necessary to find a random field without overlap, with a maximum number of iterations.
    """

    checker = True
    nr = 0

    overlap = np.zeros(2002)

    while checker:

        nr += 1
        print 'nr:',nr
        randomfieldx = randomfield(nclouds, nrbins, clon, clat, probability)

        cloud_lon_random = randomfieldx[0]
        cloud_lat_random = randomfieldx[1]
        cloud_bin_random = randomfieldx[2]

        cloud_size = binwidth * cloud_bin_random

        # check for overlap:
        cloudcentres = np.vstack((cloud_lon_random, cloud_lat_random)).T

        cloud_index = np.arange(0,nclouds)	
        pairs = itertools.combinations(cloud_index,2)

        for i in pairs:
            index1, index2 = np.array(i)
            distance = haversine(cloudcentres[index1],cloudcentres[index2])
            mindistance = cloud_size[index1]+cloud_size[index2]
            if (mindistance-distance)>0.:
                overlap[nr] = overlap[nr] + (mindistance-distance)
                #checker = True
                #break

        if overlap[nr]<binwidth/2:
            checker = False

        print overlap[nr]

        if nr>199:
            raise ValueError('Cannot find random configuration!')	

    
    return cloud_lon_random,cloud_lat_random,cloud_bin_random



def overlap(nclouds,nrbins,clon,clat,probability,binwidth):

    checker = True
    nr = 0

    overlap = np.zeros(100)
    cloud_lon_random = np.zeros([100,nclouds])
    cloud_lat_random = np.zeros([100,nclouds])
    cloud_bin_random = np.zeros([100,nclouds])

    while checker:

        #print 'nr:',nr
        randomfieldx = randomfield(nclouds, nrbins, clon, clat, probability)

        cloud_lon_random[nr,:] = randomfieldx[0]
        cloud_lat_random[nr,:] = randomfieldx[1]
        cloud_bin_random[nr,:] = randomfieldx[2]

        cloud_size = binwidth * cloud_bin_random[nr,:]

        # check for overlap:
        cloudcentres = np.vstack((cloud_lon_random[nr,:], cloud_lat_random[nr,:])).T

        cloud_index = np.arange(0,nclouds)	
        pairs = itertools.combinations(cloud_index,2)

        for i in pairs:
            index1, index2 = np.array(i)
            distance = haversine(cloudcentres[index1],cloudcentres[index2])
            mindistance = cloud_size[index1]+cloud_size[index2]
            #if (mindistance-distance)<0.:
            #    checker = False
            if (mindistance-distance)>0.:
                overlap[nr] = overlap[nr] + (mindistance - distance)       
                #checker = True
                #break

        nr += 1
        if nr>99:
            #raise ValueError('Cannot find random configuration!')	
            checker = False

   

    return cloud_lon_random,cloud_lat_random,cloud_bin_random,overlap



#randomfield = NetCDFFile('randomfield.nc','w')
#
#ncells = randomfield.createDimension('ncells',ncells_grid)
#time = randomfield.createDimension('time',0)
#height = randomfield.createDimension('height',0)
#nclouds_dim = randomfield.createDimension('nclouds_dim',nclouds)
#size_dim = randomfield.createDimension('size',0)
#
#qc = randomfield.createVariable('qc',np.float32,['time','height','ncells'])
#qc[0,0,:] = qc_fill
#qc.units = 'binary'
#qc.long_name = 'binary liquid water content (1 = water, 0 = no water'
#
#cloudlon = randomfield.createVariable('cloudlon',np.float32,['nclouds_dim'])
#cloudlon.units = 'radian'
#cloudlon.long_name = 'longitude center cloud'
#cloudlon[:] = cloud_lon_random
#
#cloudlat = randomfield.createVariable('cloud_lat',np.float64,['nclouds_dim'])
#cloudlat.units = 'radian'
#cloudlat.long_name = 'latitude center cloud'
#cloudlat[:] = cloud_lat_random
#
#cloudbin = randomfield.createVariable('cloud_bin',np.float64,['nclouds_dim'])
#cloudbin.units = '[-]'
#cloudbin.long_name = 'binnr representing size of cloud'
#cloudbin[:] = cloud_bin_random
#
#ncloudsvar = randomfield.createVariable('ncloudsvar','i4',['time'])
#ncloudsvar[:] = nclouds
#ncloudsvar.long_name = 'number of clouds in domain'
#ncloudsvar.units = '-'
#
#size_var = randomfield.createVariable('size',np.float32,['size'])
#size_var[:] = size
#
#randomfield.close()
