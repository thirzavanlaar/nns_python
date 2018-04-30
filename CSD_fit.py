#!/usr/bin/python

from matplotlib import pyplot as plt
import numpy as np
from scipy import stats 
import itertools
import math



def CSD_fit(hn,size):

    weight = np.log10(size+size[0])-np.log10(size)
 
    nb = 40
    nc = 100

    b0 = 0.
    b1 = -4.
    db = (b1-b0)/nb

    c0 = 0.
    c1 = -0.1
    dc = (c1-c0)/nc
     
    b_space = np.arange(db,nb-1,1)
    c_space = np.arange(dc,nc-1,1)    

    fitfunc = np.zeros(len(size))
    densitylog = np.where(hn>0, np.log10(hn), 0.)
    dumbias = np.zeros((nb,nc))

    combinations = list(itertools.product(b_space,c_space))
    dumbias = np.zeros((len(combinations)))
    sumsq = np.zeros((len(combinations)))

    for ib,ic in combinations:

        index = combinations.index((ib,ic))

        bfac = b0 + ib*db
        cfac = c0 + ic*dc

        fitfunc = np.exp(cfac*size)*size**bfac
        fitfunclog = np.where(fitfunc>0, np.log10(fitfunc),0.)
        nmask = np.sum(np.where((fitfunclog-densitylog)>0., 1., 0.))

        if nmask>0:
            dumbias[index] = np.sum((fitfunclog-densitylog)/nmask) #bias
            fitfunclog = fitfunclog-dumbias[index]
            sumsq[index] = np.sum(abs(weight*densitylog-fitfunclog))
            sumsq[index] = np.sum(abs(densitylog-fitfunclog))
            
    ibcmin = np.argmin(sumsq)
    
    ibmin, icmin = combinations[ibcmin]

    bfacmin = b0 + ibmin*db
    cfacmin = c0 + icmin*dc
    afacmin = dumbias[ibcmin]

    fitfunc2 = np.exp(cfacmin*size)*size**bfacmin
    logfit = np.where(fitfunc2>0., afacmin + np.log10(fitfunc2),0.)

    return afacmin, bfacmin, cfacmin, logfit
        










