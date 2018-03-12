#!/usr/bin/python

# Find minimal distances between clouds in one bin, average these per bin

from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats 
from haversine import haversine
from scipy.spatial import distance, cKDTree


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

def mindistance(filledbin,cloud_lon,cloud_lat,cloud_bin):
    ncloud_bin=np.zeros((len(size)))
    mindistance_mean=np.zeros((len(size)))
    mindistance_std=np.zeros((len(size)))

    for bb in range(0, filledbin):
        binclouds=np.zeros((nclouds,2))
        idx = np.where(cloud_bin[0,:]==bb+1)          
        binclouds[idx,0]=cloud_lon[idx]
        binclouds[idx,1]=cloud_lat[idx]
        ncloud_bin[bb] = len(idx[0])



    return mindistance_mean,mindistance_std

#mindist = np.min(spatial.distance.cdist(xy1, xy2), axis=1)
#print(mindist)


points_x = [1,2,3,4]
points_y = [2,3,4,5]

points = np.vstack((points_x,points_y)).T

Y = distance.pdist(points,'euclidean')
Z = distance.squareform(Y)

#Z = [np.nan if x==0 else x for x in Z]
Z[Z==0] = np.nan
 
print Z.shape

X = np.nanmin(Z,axis=0)
X2 = np.nanstd(Z,axis=0)

print X
print X2




        

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

