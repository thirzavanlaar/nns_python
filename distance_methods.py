#!/usr/bin/python

# Find minimal distances between clouds in one bin, average these per bin

from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats 
from haversine import haversine
from scipy.spatial import distance


cusize = NetCDFFile(
         '/home/vanlaar/HDCP2data/TA_dom4/cusize_output_time41.nc')

cloudlon = cusize.variables['cloud_lon']
cloudlat = cusize.variables['cloud_lat']
nclouds_cusize  = cusize.variables['nclouds']
size = cusize.variables['size']
cloud_bin = cusize.variables['cloud_bin']
hn = cusize.variables['hn']

nclouds = int(nclouds_cusize[0])
cloud_lon = cloudlon[0,0:nclouds]
cloud_lat = cloudlat[0,0:nclouds]
filledbin=np.argmin(hn[0,:])

def binclouds(cloud_lon,cloud_lat,filledbin,nclouds):
    ncloud_bin=np.zeros((len(size)))
    binclouds=np.zeros((nclouds,2))
    for bb in range(0, filledbin):
        idx = np.where(cloud_bin[0,:]==bb+1)          
        binclouds[idx,0]=cloud_lon[idx]
        binclouds[idx,1]=cloud_lat[idx]
        ncloud_bin[bb] = len(idx[0])
    return [binclouds,ncloud_bin]


def binclouds_old(cloud_lon,cloud_lat,filledbin,nclouds):
    
    ncloud_bin=np.zeros((len(size)))
    binclouds=np.zeros((nclouds,2))

    for bb in range(0, filledbin):
    #for bb in range(6, 7):
        print 'bin',bb
        cc=0
        for aa in range(1,nclouds):
            if cloud_bin[0,aa] == bb+1:
                binclouds[cc,0]=cloud_lon[aa]
                binclouds[cc,1]=cloud_lat[aa]
                cc=cc+1
        ncloud_bin[bb] = int(cc)

    return [binclouds,ncloud_bin]

def mindistance_old(binclouds,filledbin):
    
    dist=np.zeros((nclouds*(nclouds-1)/2))
    mindistance_mean=np.zeros((len(size)))
    mindistance_std=np.zeros((len(size)))
    
    for bb in range(0, filledbin):
    #for bb in range(6, 7):
        cc = int(ncloud_bin[bb])
        for ff in range(0, cc):
            distold = 400000.
            disttotal = 0.
            firstx = binclouds[ff,0]
            firsty = binclouds[ff,1]
            for ss in range(0,cc):
                if ss!=ff:
                    secondx = binclouds[ss,0]
                    secondy = binclouds[ss,1]
                    distnew = haversine(firsty,firstx,secondy,secondx)
                    if distnew<=distold:
                        distold = distnew
            dist[ff] = distold
        mean_distances = np.mean(dist[0:cc])
        std_distances = np.std(dist[0:cc])
        if cc>2:
            mindistance_mean[bb] = mean_distances
        else:
            mindistance_mean[bb] = 0.
        if cc>3:
            mindistance_std[bb] = std_distances
        else:
            mindistance_std[bb] = 0

    return mindistance_mean,mindistance_std      


def meandistance(filledbin,cloud_lon,cloud_lat,cloud_bin):
    ncloud_bin=np.zeros((len(size)))
    D1=np.zeros((len(size)))
    D0=np.zeros((len(size)))

    for bb in range(0, filledbin):
        binclouds=np.zeros((nclouds,2))
        idx = np.where(cloud_bin[0,:]==bb+1)          
        binclouds[idx,0]=cloud_lon[idx]
        binclouds[idx,1]=cloud_lat[idx]
        ncloud_bin[bb] = len(idx[0])
        Y = distance.pdist(binclouds[idx],haversine)
        D0[bb] = stats.gmean(Y)
        D1[bb] = np.mean(Y)
    return D0,D1
        
output_meandistance = meandistance(filledbin,cloud_lon,cloud_lat,cloud_bin)
D0 = output_meandistance[0]
D1 = output_meandistance[1]


#output_mindistance = mindistance(binclouds,filledbin)
#mindistance_mean = output_mindistance[0]
#mindistance_std = output_mindistance[1]
#
#
#mindistance_total = mindistance_mean+mindistance_std
#
#plt.figure(figsize=(10,8))
#plt.axis([0, 6000, 0, 70000])
#plt.xlabel('Size')
#plt.ylabel('Minimal distance')
#plt.scatter(size,mindistance_mean,color='k')
#plt.scatter(size,mindistance_total,color='g')
#plt.savefig('mindistance.pdf')


plt.figure(figsize=(10,8))
#plt.axis([50000,220000, 50000, 220000])
plt.xlabel('D1')
plt.ylabel('D0')
plt.scatter(D1,D0,color='k')
plt.savefig('Figures/D1-D0.pdf')

