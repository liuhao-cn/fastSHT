#!/usr/bin/env python
# coding: utf-8

# In[1]:
import sys as sys
import os

# the command line input will overwrite the defaults
nside = int(sys.argv[1])
nsim = int(sys.argv[2])
n_proc = int(sys.argv[3])
niter = int(sys.argv[4])


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

import importlib
importlib.reload(SHT)


nrep = 1
lmax = 3*nside - 1
npix = 12 * nside ** 2

def test_t2alm(seed=23333):

    np.random.seed(seed)
    # maps = np.asfortranarray(np.random.rand(npix, nsim))
    maps = np.ones([npix, nsim], dtype=np.double, order='F')
    alms = np.ones([nsim, lmax+1, lmax+1], dtype=np.double,  order='F')

    sht = SHT.SHT(nside, lmax, nsim, niter)

    start = time.time()
    sht.t2alm(maps, alms_in=alms)
    end1 = time.time() - start

    start = time.time()
    for i in range(nsim):
        hp.map2alm(maps[:,i], lmax=lmax, iter=niter)
    end2 = time.time() - start
    
    return end1, end2

end1, end2 = test_t2alm()


filename = '%4.4i-%5.5i.txt' %(nside, nsim)
with open(filename, 'a') as f:
    sys.stdout = f
    print("t2alm: Nside = %i, Nsim = %i, n_proc = %2i, fastSHT-CPU = %6.2f, healpy = %6.2f" 
        %(nside, nsim, n_proc, end1, end2))

