#!/usr/bin/env python
# coding: utf-8

# Import modules

import sys as sys
import os

seed = 23333

n_proc = 8
os.environ["OMP_NUM_THREADS"] = str(n_proc) # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = str(n_proc) # export OPENBLAS_NUM_THREADS=4
os.environ["MKL_NUM_THREADS"] = str(n_proc) # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = str(n_proc) # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = str(n_proc) # export NUMEXPR_NUM_THREADS=6

sys.path.append('..')
from SHT import SHT

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt



# Set main parameters

nside = 128        # nside
lmax = 3*nside-1   # lmax
nsim = 100         # number of maps
niter = 3          # number of iterations

npix = 12*nside**2



# Creating test temperature and polarization maps 

np.random.seed(seed)

T = np.asfortranarray(np.random.rand(npix, nsim))
Q = np.asfortranarray(np.random.rand(npix, nsim))
U = np.asfortranarray(np.random.rand(npix, nsim))



# Initialize the SHT handle

# For t2alm, both pol=False and pol=True works fine, but pol=False will use less memory
sht = SHT(nside, lmax, nsim, niter, pol=False)



# T to alms by fastSHT

alms = sht.t2alm(T)



# T to alms by healpy

# convert alms into the healpy format
alms_hp1 = sht.convert_alm_healpy(alms)

# compute alms directly from Healpy
alms_hp2 = alms_hp1 * 0

for i in range(nsim):
    alms_hp2[:,i] = hp.map2alm(T[:,i], lmax=lmax, iter=niter)



# Test the TT-accuracy with Healpy

cl1 = np.array([hp.alm2cl(alms_hp1[:,i]) for i in range(nsim)])
cl2 = np.array([hp.alm2cl(alms_hp2[:,i]) for i in range(nsim)])

print('Max. relative TT-error:', (np.abs(cl2-cl1)/cl2.mean()).max() )



# QU to EB

# For QU2EB one has to use pol=True
sht = SHT(nside, lmax, nsim, niter, pol=True)

almEs, almBs = sht.qu2eb(Q, U)

almEs_hp1 = sht.convert_alm_healpy(almEs)
almBs_hp1 = sht.convert_alm_healpy(almBs)

# test with Healpy for EE and BB accuracy
maps = np.asfortranarray( [T, Q, U ] )

alms_hp2 = np.array([hp.map2alm(maps[:,:,i], lmax=lmax, iter=niter) for i in range(nsim)])
almEs_hp2 = alms_hp2[:,1,:].transpose()
almBs_hp2 = alms_hp2[:,2,:].transpose()

cl1 = np.array([hp.alm2cl(almEs_hp1[:,i]) for i in range(nsim)])
cl2 = np.array([hp.alm2cl(almEs_hp2[:,i]) for i in range(nsim)])
max_errEE = (np.abs(cl2 - cl1) / cl2.mean()).max()
print('Max. relative EE-error:', max_errEE )

cl1 = np.array([hp.alm2cl(almBs_hp1[:,i]) for i in range(nsim)])
cl2 = np.array([hp.alm2cl(almBs_hp2[:,i]) for i in range(nsim)])
max_errBB = (np.abs(cl2 - cl1) / cl2.mean()).max()
print('Max. relative BB-error:', max_errBB )



# Alms to T

# convert alms to maps
T1 = sht.alm2t(alms)

T2 = np.array([hp.alm2map(alms_hp2[i,0,:], nside, pol=False) for i in range(nsim)]).transpose()

max_errT = (np.abs(T2 - T1)/T2.std()).max()

print('Max. relative T-error:', max_errT )



# EB to QU

Q1, U1 = sht.eb2qu(almEs, almBs)

maps2 = np.array([hp.alm2map(alms_hp2[i,:,:], nside, pol=True) for i in range(nsim)])
Q2 = maps2[:,1,:].transpose()
U2 = maps2[:,2,:].transpose()

max_errQ = (np.abs(Q2 - Q1)/Q1.std()).max()
max_errU = (np.abs(U2 - U1)/U1.std()).max()

print('Max. relative Q-error:', max_errQ )
print('Max. relative U-error:', max_errU )


