#!/usr/bin/env python
# coding: utf-8

# In[1]:
import sys as sys
import os

nside = 64
nsim = 2000
n_proc = 8
niter = 3
compare = False

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
print("nside = %i, nsim = %i, n_proc = %i, niter = %i, comparison = %s" 
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

import numba
from numba import cuda


# In[2]:


import importlib
importlib.reload(SHT)

# In[3]:


lmax = 3*nside - 1
npix = 12 * nside ** 2
test_cl = np.array([1 for l in range(1,lmax+1)] )

def test_t2alm(nside, lmax, nsim, test_cl, niter = 0, seed=23333, compare=True):
    print('Testing t2alm...')

    np.random.seed(seed)
    maps = np.asfortranarray(np.transpose([hp.sphtfunc.synfast(test_cl, nside, lmax)
                 for i in range(nsim)]) )

    #maps = np.load('maps_' + str(nsim) + '.npy')

    start = time.time()
    sht = SHT.SHT(nside, lmax, nsim, niter)
    #alms = numba.cuda.pinned_array((nsim, lmax+1, lmax+1), dtype=np.double, strides=None, order='F')
    alms = np.empty((nsim, lmax+1, lmax+1), dtype=np.double,  order='F')
    end = time.time() - start
    print('Time cost for memory initialization is ' + str(end))

    start = time.time()
    for i in range(5):
        sht.t2alm_old(maps, alms)
    end = time.time() - start
    print('Calculation time cost for fastSHT is ' + str(end / 5))

    #print(time.sleep(10))
    if(compare==False):
        return
    alms_hp = sht.convert_alm_healpy(alms)
    alms_hp = (alms_hp[0,:,:] + 1j * alms_hp[1,:,:])

    del alms
    alm_shape = alms_hp.shape[0]
    del alms_hp
    
    start = time.time()
    overhead = 0
    roverhead = 0
    #alms2_hp = np.array([hp.map2alm(maps[:,i], lmax=lmax, iter=niter) for i in range(nsim)])
    for j in range(5):
        for i in range(nsim):
            ostart = time.time()
            maps[:,i].astype(np.float64, order='C', copy=True)
            rstart = time.time()
            np.empty((alm_shape),order='C', dtype=np.complex128)
            end = time.time()
            overhead += end - ostart
            roverhead += end - rstart
            hp.map2alm(maps[:,i], lmax=lmax, iter=niter)

    tot_time = time.time() - start
    print('Time cost for memory initialization is ' + str(roverhead/5))
    print('Calculation time cost for healpy is ' + str((tot_time - overhead)/5 ))
    #print('Total time cost for healpy is ' + str(tot_time/5))
    
    # cl = np.array([hp.alm2cl(alms_hp[:,i]) for i in range(nsim)])
    # cl2 = np.array([hp.alm2cl(alms2_hp[i,:]) for i in range(nsim)])
    
    # max_err = (np.abs(cl2 - cl) / cl2.mean()).max()

    # print('Max relative error in cl is: ' + str(max_err))


def make_mask(nside, upper_lat = 5):
    npix = nside**2 * 12
    mask = np.ones(npix)
    
    for i in range(npix):
        theta, phi = hp.pix2ang(nside, i, lonlat=True)
        if phi < 5 and phi > -5:
            mask[i] = 0
    return mask

# In[7]:

def fix_EB(Q, U, mask, nside, lmax, niter=0, seed=23333):
    
            
    vid = (np.arange(len(mask))[mask == 1])
    nv = len(vid)

    maps_in = np.array( [np.transpose([hp.sphtfunc.synfast(test_cl, nside, lmax) for i in range(nsim)]), Q, U ])
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
        
        # mx = np.sum(x) / nv
        # my = np.sum(y) / nv
        # cxx = np.sum( (x - mx)**2 )
        # cxy = np.sum( (y - my) * (x - mx) )
        # coe = [cxy / cxx, my  - mx * (cxy / cxx) ]
        
        BC[i, vid] = BO[i, vid] - BT[i, vid] * coe[0] - coe[1]
        
    print('Time cost for Healpy is ' + str(time.time() - start))
    return BC

def test_fix_EB(nside, lmax, niter=0, seed=23333, compare=True):
    npix = 12*nside**2
    print('Testing fix_EB...')
    np.random.seed(seed)
    Q = np.asfortranarray(np.transpose([hp.sphtfunc.synfast(test_cl, nside, lmax)
                 for i in range(nsim)]) )
    U = np.asfortranarray(np.transpose([hp.sphtfunc.synfast(test_cl, nside, lmax)
                 for i in range(nsim)]) )
    
    mask = make_mask(nside)

    start = time.time()
    sht = SHT.SHT(nside, lmax, nsim, niter, pol=True)

    Bmap = sht.fix_eb(Q, U, mask)
    print('Time cost for fastSHT is ' + str(time.time() - start))
    if(compare == False):
        return
    
    
    BC = fix_EB(Q, U, mask, nside, lmax, niter, seed)

    cl = np.array([hp.anafast(Bmap[:,i]) for i in range(nsim)])
    cl2 = np.array([hp.anafast(BC[i,:]) for i in range(nsim)])
    max_err = (np.abs(cl2 - cl) / cl.mean()).max()
    print('Max relative error in clB is: ' + str(max_err))    
    return (Bmap, BC)

def test_qu2eb(nside, lmax, nsim, test_cl, niter = 0, seed=23333, compare=True):
    print('Testing qu2eb...')
    np.random.seed(seed)
    Q = np.asfortranarray(np.transpose([hp.sphtfunc.synfast(test_cl, nside, lmax)
                 for i in range(nsim)]) )
    U = np.asfortranarray(np.transpose([hp.sphtfunc.synfast(test_cl, nside, lmax)
                 for i in range(nsim)]) )
    
    
    sht = SHT.SHT(nside, lmax, nsim, niter, pol=True)
    
    start = time.time()
    almEs, almBs = sht.qu2eb(Q, U)
    print('Time cost for fastSHT is ' + str(time.time() - start))

    if(compare == False):
        return

    almEs_hp = sht.convert_alm_healpy(almEs)
    almEs_hp = (almEs_hp[0,:,:] + 1j * almEs_hp[1,:,:])
    almBs_hp = sht.convert_alm_healpy(almBs)
    almBs_hp = (almBs_hp[0,:,:] + 1j * almBs_hp[1,:,:])
    
    
    maps = np.asfortranarray( [np.transpose([hp.sphtfunc.synfast(test_cl, nside, lmax) for i in range(nsim)]), Q, U ])
    
    start = time.time()
    alms2_hp = np.array([hp.map2alm(maps[:,:,i], lmax=lmax, iter=niter) for i in range(nsim)])
    print('Time cost for healpy is ' + str(time.time() - start))
    
    cl = np.array([hp.alm2cl(almEs_hp[:,i]) for i in range(nsim)])
    cl2 = np.array([hp.alm2cl(alms2_hp[i, 1, :]) for i in range(nsim)])
    
    max_err = (np.abs(cl2 - cl) / cl.mean()).max()
    print('Max relative error in clE is: ' + str(max_err))
    
    cl = np.array([hp.alm2cl(almBs_hp[:,i]) for i in range(nsim)])
    cl2 = np.array([hp.alm2cl(alms2_hp[i, 2, :]) for i in range(nsim)])
    
    max_err = (np.abs(cl2 - cl) / cl.mean()).max()
    print('Max relative error in clB is: ' + str(max_err))
# In[6]:

    
#test_t2alm(nside, lmax, nsim, np.array([1 for l in range(1,lmax+1)] ), niter=niter, compare=compare)

#test_qu2eb(nside, lmax, nsim, np.array([1 for l in range(1,lmax+1)] ), niter=niter, compare=compare)
test_fix_EB(nside, lmax, niter, compare=compare)
