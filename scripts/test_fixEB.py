#!/usr/bin/env python
# coding: utf-8

# In[1]:
import sys as sys
import os

nside = 64
nsim = 1000
n_proc = 8
niter = 3
compare = False

# the command line input will overwrite the defaults
if len(sys.argv)>1:
    nside = int(sys.argv[1])
if len(sys.argv)>2:
    nsim = int(sys.argv[2])
if len(sys.argv)>3:
    n_proc = int(sys.argv[3])
if len(sys.argv)>4:
    niter = int(sys.argv[4])
if len(sys.argv)>5:
    compare = sys.argv[5].lower() == "true"

print(" ")
print("Working with the following parameters:")
print("Nside = %i, Nsim = %i, n_proc = %i, Niter = %i, comparison with HEALPix = %s" 
    %(nside, nsim, n_proc, niter, compare))
print(" ")

os.environ["OMP_NUM_THREADS"] = str(n_proc) # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = str(n_proc) # export OPENBLAS_NUM_THREADS=4 
os.environ["MKL_NUM_THREADS"] = str(n_proc) # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = str(n_proc) # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = str(n_proc) # export NUMEXPR_NUM_THREADS=6

sys.path.append('..')
import SHT
import fastSHT


import numpy as np

import healpy as hp
import matplotlib.pyplot as plt
import time

import importlib
importlib.reload(SHT)


def make_mask(nside, upper_lat = 5):
    npix = nside**2 * 12
    mask = np.ones(npix)
    
    for i in range(npix):
        theta, phi = hp.pix2ang(nside, i, lonlat=True)
        if phi < 5 and phi > -5:
            mask[i] = 0
    return mask

def fix_EB_hp(T, Q, U, mask, seed=23333):
    
    vid = (np.arange(len(mask))[mask == 1])
    nv = len(vid)

    maps_in = np.array( [T, Q, U ] )
    maps_in[:,mask==0,:] = hp.UNSEEN
    
    start = time.time()

    alms2_hp = np.array([hp.map2alm(maps_in[:,:,i], lmax=lmax, iter=niter) for i in range(nsim)])
    
    BO = np.array([hp.alm2map(alms2_hp[i,2,:], nside=nside, lmax=lmax) for i in range(nsim)])
    
    alms2_hp[:,2,:] = 0
    alms2_hp[:,0,:] = 0
    
    maps = np.array([hp.alm2map(alms2_hp[i,:,:], nside=nside, lmax=lmax) for i in range(nsim)])
    maps[:,:,mask==0] = 0
    alms2_hp = np.array([hp.map2alm(maps[i,:,:], lmax=lmax, iter=niter) for i in range(nsim)])
    alms2_hp = alms2_hp[:,2,:]
    BT = np.array([hp.alm2map(alms2_hp[i,:], nside=nside, lmax=lmax) for i in range(nsim)])
    BC = np.zeros(BT.shape)
    for i in range(nsim):
        x = BT[i, vid]
        y = BO[i, vid]
        coe = np.polyfit(x, y, 1)

        BC[i, vid] = BO[i, vid] - BT[i, vid] * coe[0] - coe[1]
        
    print('Time cost for Healpy is ' + str(time.time() - start))
    return BC
    


def test_fix_EB(seed=23333):

    print('Testing fix_EB...')
    np.random.seed(seed)

    T = np.asfortranarray(np.random.rand(npix, nsim))
    Q = np.asfortranarray(np.random.rand(npix, nsim))
    U = np.asfortranarray(np.random.rand(npix, nsim))
    
    mask = make_mask(nside)
    
    sht = SHT.SHT(nside, lmax, nsim, niter, pol=True)

    start = time.time()
    Bmap = sht.fix_eb(Q, U, mask)
    print('Time cost for fastSHT is ' + str(time.time() - start))
    if(compare == False):
        return
    
    BC = fix_EB_hp(T, Q, U, mask, seed=seed)

    cl = np.array([hp.anafast(Bmap[:,i]) for i in range(nsim)])
    cl2 = np.array([hp.anafast(BC[i,:]) for i in range(nsim)])
    max_err = (np.abs(cl2 - cl) / cl.mean()).max()
    print('Max relative error in clB is: ' + str(max_err))    
    return (Bmap, BC)
        


npix = 12*nside**2
lmax = 3*nside - 1

test_fix_EB()


