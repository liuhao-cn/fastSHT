#!/usr/bin/env python
# coding: utf-8

# In[1]:
import sys as sys
import os
n_proc = 8
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


# In[2]:


import importlib
importlib.reload(SHT)

# In[3]:


nside = 128
lmax = 383
nsim = 200
npix = 12 * nside ** 2
test_cl = np.array([1 for l in range(1,lmax+1)] )


# In[ ]:




# In[4]:


def test_t2alm(nside, lmax, nsim, test_cl, niter = 0, seed=23333, compare=True):
    print('Testing t2alm...')
    np.random.seed(seed)
    maps = np.asfortranarray(np.transpose([hp.sphtfunc.synfast(test_cl, nside, lmax)
                 for i in range(nsim)]) )
    #maps = np.save('maps_' + str(nsim) + '.npy', maps)
    #maps = np.load('maps_' + str(nsim) + '.npy')
    sht = SHT.SHT(nside, lmax, nsim, niter)
    
    start = time.time()
    alms = sht.t2alm_old(maps)
    print('Time cost for fastSHT is ' + str(time.time() - start))

    if(compare == False):
        return
    
    alms_hp = sht.convert_alm_healpy(alms)
    alms_hp = (alms_hp[0,:,:] + 1j * alms_hp[1,:,:])
    
    start = time.time()
    alms2_hp = np.array([hp.map2alm(maps[:,i], lmax=lmax, iter=niter) for i in range(nsim)])
    print('Time cost for healpy is ' + str(time.time() - start))
    
    cl = np.array([hp.alm2cl(alms_hp[:,i]) for i in range(nsim)])
    cl2 = np.array([hp.alm2cl(alms2_hp[i,:]) for i in range(nsim)])
    
    max_err = (np.abs(cl2 - cl) / cl2.mean()).max()

    #print(cl[0,0:10])
    #print(cl2[0,0:10])
    print('Max relative error in cl is: ' + str(max_err))
# In[5]:


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
    
    sht = SHT.SHT(nside, lmax, nsim, niter, pol=True)
    start = time.time()
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
        

# In[ ]:




# In[9]:


nsim = 200
test_t2alm(nside, lmax, nsim, np.array([1 for l in range(1,lmax+1)] ), niter=0, compare=True)

test_t2alm(nside, lmax, nsim, np.array([1 for l in range(1,lmax+1)] ), niter=3, compare=True)

test_qu2eb(nside, lmax, nsim, np.array([1 for l in range(1,lmax+1)] ), niter=0, compare=True)

test_qu2eb(nside, lmax, nsim, np.array([1 for l in range(1,lmax+1)] ), niter=3, compare=True)

nside = 64
lmax = 191
nsim = 200
npix = 12 * nside ** 2
test_cl = np.array([1 for l in range(1,lmax+1)] )

test_fix_EB(nside, lmax, 3, compare=True)
nsim = 2000
#test_t2alm(nside, lmax, nsim, np.array([1 for l in range(1,lmax+1)] ), niter=0, compare=True)

#test_t2alm(nside, lmax, nsim, np.array([1 for l in range(1,lmax+1)] ), niter=3, compare=True)

#test_qu2eb(nside, lmax, nsim, np.array([1 for l in range(1,lmax+1)] ), niter=0, compare=False)

# In[5]:


