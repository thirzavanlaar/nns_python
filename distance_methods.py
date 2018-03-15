#!/usr/bin/python

# Find minimal distances between clouds in one bin, average these per bin
# Compute geometric and arithmetical mean between all clouds per bin

from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats 
from haversine import haversine
from scipy.spatial import distance

# Import data:

#cusize = NetCDFFile(
#         '/home/vanlaar/HDCP2data/TA_dom4/cusize_output_time41.nc')
#
#cloudlon = cusize.variables['cloud_lon']
#cloudlat = cusize.variables['cloud_lat']
#nclouds_cusize  = cusize.variables['nclouds']
#size = cusize.variables['size']
#cloud_bin = cusize.variables['cloud_bin']
#hn = cusize.variables['hn']
#
#nclouds = int(nclouds_cusize[0])
#cloud_lon = cloudlon[0,0:nclouds]
#cloud_lat = cloudlat[0,0:nclouds]
#filledbin=np.argmin(hn[0,:])  # last bin with clouds, rest is empty


def distances(filledbin,cloud_lon,cloud_lat,cloud_bin,size,nclouds):
    ncloud_bin=np.zeros((len(size)))
    D1=np.zeros((len(size)))
    D0=np.zeros((len(size)))
    mindistance_mean=np.zeros((len(size)))
    mindistance_std=np.zeros((len(size)))

    for bb in range(0, filledbin+1):
        binclouds=np.zeros((nclouds,2))
        # Select clouds present in current bin bb:
        idx = np.where(cloud_bin[0,:]==bb+1)          
        binclouds[idx,0]=cloud_lon[idx]
        binclouds[idx,1]=cloud_lat[idx]
        ncloud_bin[bb] = len(idx[0])
        if ncloud_bin[bb]>=3.:
            # Compute all distances by using the haversine method:
            Y = distance.pdist(binclouds[idx],haversine)
            D0[bb] = stats.gmean(Y)
            D1[bb] = np.mean(Y)

            # Taking only the minimum distance per cloud (nearest neighbours):
            Z = distance.squareform(Y)
            Z[Z==0] = np.nan
            mindistances = np.nanmin(Z,axis=0)
            mindistance_mean[bb] = np.mean(mindistances)
            mindistance_std[bb] = np.std(mindistances)
    return D0,D1,mindistance_mean,mindistance_std,ncloud_bin
        

def JosephCahalan(filledbin,cloud_bin,cloud_lon,cloud_lat,size):

    cloudcentres = np.vstack((cloud_lon,cloud_lat)).T 
    
    mindistance_JC_mean = np.zeros((len(size)))
    mindistance_JC_std = np.zeros((len(size)))

    Y = distance.pdist(cloudcentres,haversine)
    Z = distance.squareform(Y)
    Z[Z==0] =np.nan
    mindistances = np.nanmin(Z,axis=0)

    for bb in range(0, filledbin+1):
        idx = np.where(cloud_bin[0,:]==bb+1)
        if len(idx[0]) >= 3:
            mindistance_JC_mean[bb] = np.nanmean(mindistances[idx])
            mindistance_JC_std[bb] = np.nanstd(mindistances[idx])

    return mindistance_JC_mean,mindistance_JC_std













