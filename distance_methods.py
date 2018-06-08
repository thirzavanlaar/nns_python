#!/usr/bin/python

# Find minimal distances between clouds in one bin, average these per bin
# Compute geometric and arithmetical mean between all clouds per bin

from netCDF4 import Dataset as NetCDFFile
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats 
from haversine import haversine, haversine_V
from scipy.spatial import distance
from collections import Counter



def distances(filledbin,cloud_lon,cloud_lat,cloud_bin,size,nclouds):
    """
    Compute nearest neighbour distances only between clouds belonging to the same bin.
    """
    ncloud_bin=np.zeros((len(size)))
    D1=np.zeros((len(size)))
    D0=np.zeros((len(size)))
    mindistance_mean=np.zeros((len(size)))
    mindistance_std=np.zeros((len(size)))
    maxdistance=np.zeros((len(size)))
    maxdistanceY=np.zeros((len(size)))

    for bb in range(0, filledbin+1):
        binclouds=np.zeros((nclouds,2))
        # Select clouds present in current bin bb:
        idx = np.where(cloud_bin[:]==bb+1)          
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
            Z = np.ma.masked_where(Z==0,Z)
            mindistances = np.min(Z,axis=0)
            mindistance_mean[bb] = np.ma.mean(mindistances)
            mindistance_std[bb] = np.std(mindistances)
            maxdistance[bb] = np.max(mindistances)
            maxdistanceY[bb] = np.max(Y)
    return D0,D1,mindistance_mean,mindistance_std,ncloud_bin
        

def JosephCahalan(filledbin,cloud_bin,cloud_lon,cloud_lat,size):
    """
    Compute nearest neighbour distances between all clouds and average them per bin, based on Joseph and Cahalan (1990).
    """

    cloudcentres = np.vstack((cloud_lon,cloud_lat)).T 

    nr_bins = len(size)    
    ncloud_bin=np.zeros((len(size)))
    mindistance_JC_mean = np.zeros((len(size)))
    mindistance_JC_std = np.zeros((len(size)))
    neighbour_avg = np.zeros((len(size)))
    neighbour_histo = np.zeros((len(size),len(size)))

    Y = distance.pdist(cloudcentres,haversine)
    Z = distance.squareform(Y)
    Z = np.ma.masked_where(Z==0,Z)
    mindistances = np.min(Z,axis=0)

    neighbour_idx = np.argmin(Z,axis=0)
    neighbour_size = np.array([cloud_bin[0,xi] for xi in neighbour_idx])

    bins_histo = np.arange(0,len(size)+1)

    for bb in range(0, filledbin+1):
        idx = np.where(cloud_bin[0,:]==bb+1)
        ncloud_bin[bb] = len(idx[0])
        if len(idx[0]) >= 3:
            mindistance_JC_mean[bb] = np.mean(mindistances[idx])
            mindistance_JC_std[bb] = np.std(mindistances[idx])          
            neighbour_avg[bb] = np.mean(neighbour_size[idx])
            histogram = np.histogram(neighbour_size[idx],bins=bins_histo,density=True)
            neighbour_histo[bb,:] = histogram[0]


    return mindistance_JC_mean,mindistance_JC_std,neighbour_avg,neighbour_histo,ncloud_bin


def JosephCahalan_kneighbour(filledbin,cloud_bin,cloud_lon,cloud_lat,size,nr_neighbours):
    """
    Take the k nearest neighbours of a cloud, average the k distances to these clouds and then average these values per bin.
    """

    cloudcentres = np.vstack((cloud_lon,cloud_lat)).T 

    nr_bins = len(size)    
    mindistance_JC_mean = np.zeros((len(size)))
    mindistance_JC_std = np.zeros((len(size)))
    neighbour_avg = np.zeros((len(size)))
    neighbour_histo = np.zeros((len(size),len(size)))

    Y = distance.pdist(cloudcentres,haversine)
    Z = distance.squareform(Y)
    Z = np.ma.masked_where(Z==0,Z)

    mindistances = np.zeros(len(cloud_lon))
    neighbour_size = np.zeros(len(cloud_lon))

    for i in range(0,len(cloud_lon)):
        A = Z[i,:]
        idx = np.argpartition(A,nr_neighbours)
        kneighbours = A[idx[:nr_neighbours]]
        mindistances[i] = np.mean(kneighbours)
        ksizes = cloud_bin[0,idx[:nr_neighbours]]
        neighbour_size[i] = np.mean(ksizes)

    #neighbour_idx = np.argpartition(

    #mindistances = np.min(Z,axis=0)

    #neighbour_idx = np.argmin(Z,axis=0)
    #neighbour_size = np.array([cloud_bin[0,xi] for xi in neighbour_idx])

    #histo = Counter(neighbour_size)
    #histo = np.histogram(neighbour_size)
    
    bins_histo = np.arange(0,len(size)+1)

    for bb in range(0, filledbin+1):
        idx = np.where(cloud_bin[0,:]==bb+1)
        #print 'idx:',idx
        if len(idx[0]) >= 3:
            mindistance_JC_mean[bb] = np.mean(mindistances[idx])
            mindistance_JC_std[bb] = np.std(mindistances[idx])          
            neighbour_avg[bb] = np.mean(neighbour_size[idx])
            histogram = np.histogram(neighbour_size[idx],bins=bins_histo,density=True)
            #len_hist = len(histogram[0])
            neighbour_histo[bb,:] = histogram[0]


    return mindistance_JC_mean,mindistance_JC_std,neighbour_avg,neighbour_histo

def SCAI(cloud_lon,cloud_lat,N,Nmax,L):

    """
    Compute the SCAI for a cloud field, based on Tobin et al (2012). 
    """

    cloudcentres = np.vstack((cloud_lon,cloud_lat)).T 

    Y = distance.pdist(cloudcentres,haversine)
    
    D0 = stats.gmean(Y)
    D1 = np.mean(Y)
   
    N = float(N)
    Nmax = float(Nmax)
 
    SCAI_0 = (N/Nmax)*(D0/L)*1000
    SCAI_1 = (N/Nmax)*(D1/L)*1000

    return SCAI_0,SCAI_1



def COP(cloud_lon,cloud_lat,cloud_size):
    """
    Compute the COP for a cloud field, based on White et al (2017).
    """

    cloudcentres = np.vstack((cloud_lon,cloud_lat,cloud_size)).T 

    V = distance.pdist(cloudcentres,haversine_V)
    
    COP = np.sum(V)/len(V)

    return COP


def COP_perbin(cloud_lon,cloud_lat,cloud_bin,size,filledbin):
    """
    Compute the COP for a cloud field, based on White et al (2017).
    """

    cloud_size = cloud_bin*size[0]
    cloudcentres = np.vstack((cloud_lon,cloud_lat,cloud_size)).T 

    COP_perbin = np.zeros(len(size))

    for bb in range(0, filledbin+1):
        idx = np.where(cloud_bin[:]==bb+1)
        cloudsinbin = cloudcentres[idx,:]
        V = distance.pdist(cloudsinbin[0,:,:],haversine_V)
        COP_perbin[bb] = np.sum(V)/len(V)
        #COP_perbin[bb] = np.sum(V)
        
    return COP_perbin 



def NNCDF(cloud_lon,cloud_lat):
    """
    Compute the nncdf of a cloud field, based on Nair et al (1998).
    """

    cloudcentres = np.vstack((cloud_lon,cloud_lat)).T 

    Y = distance.pdist(cloudcentres,haversine)
    Z = distance.squareform(Y)
    Z = np.ma.masked_where(Z==0,Z)
    mindistances = np.min(Z,axis=0)

    distance_max = max(mindistances)
    distance_min = min(mindistances)    

    print distance_max/100.
    print distance_max

    step = distance_max/100.

    print mindistances

    bin_edges = np.arange(0, distance_max+step, step)
    print bin_edges.shape
    #bins_histo = np.arange(0,filledbin)
    #bins_histo = np.arange(0,len(size)+1)
    #bins_histo = np.arange(0,100)
    values, base = np.histogram(mindistances,bins=bin_edges,density=True)
    print sum(values)

    nncdf = np.cumsum(values)

    print nncdf
    return nncdf, base


def distance_poisson(filledbin,cloud_bin,cloud_lon,cloud_lat,size):
    """
    """

    cloudcentres = np.vstack((cloud_lon,cloud_lat)).T 

    nr_bins = len(size)    
    ncloud_bin=np.zeros((len(size)))
    mindistance_JC_mean = np.zeros((len(size)))
    mindistance_JC_std = np.zeros((len(size)))
    neighbour_avg = np.zeros((len(size)))
    neighbour_histo = np.zeros((len(size),len(size)))

    Y = distance.pdist(cloudcentres,haversine)
    Z = distance.squareform(Y)
    Z = np.ma.masked_where(Z==0,Z)
    mindistances = np.min(Z,axis=0)

    print mindistances.shape

    return mindistances,ncloud_bin


