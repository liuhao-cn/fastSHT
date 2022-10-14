
import sys
from pathlib import Path
path = Path(__file__).parent 
sys.path.append(str(path.resolve()) + '/obj/')
import fastSHT

import numpy as np
import healpy as hp

import time

import numba
from numba import cuda

class SHT:
    def __init__(self, nside, lmax, nsim, niter=0, pol=False, all_buff=False):
        self.nside = nside
        self.lmax = lmax
        self.nsim = nsim
        self.niter = niter
        if(pol == True):
            self.pol = 1
        else:
            self.pol = 0


        self.nring = 4 * nside - 1
        self.npix = 12 * nside ** 2
        
        ring_tab_dir = str(path.resolve()) + '/SHT_ring_skip_table/'
        ring_tab = np.fromfile(ring_tab_dir + str(nside).rjust(4,'0') + '.bin', dtype=np.int32)

        plm_dir = str(path.resolve()) + '/plm_table/'
        plm_pos = np.fromfile(plm_dir + 'plm_pos_' + str(nside).rjust(4,'0') + '.bin', dtype=np.int32)
        plm_val1 = np.fromfile(plm_dir + 'plm_val1_' + str(nside).rjust(4,'0') + '.bin', dtype=np.double)
        plm_val2 = np.fromfile(plm_dir + 'plm_val2_' + str(nside).rjust(4,'0') + '.bin', dtype=np.double)

        plm_pos = np.asfortranarray(plm_pos.reshape(2*nside, int(len(plm_pos) / (2*nside)), order='F')[:,:lmax+1], dtype=np.int32)
        plm_val1 = np.asfortranarray(plm_val1.reshape(2*nside, int(len(plm_val1) / (2*nside)), order='F')[:,:lmax+1], dtype=np.double)
        plm_val2 = np.asfortranarray(plm_val2.reshape(2*nside, int(len(plm_val2) / (2*nside)), order='F')[:,:lmax+1], dtype=np.double)

        plm_val1[plm_pos >= lmax] = 0
        plm_val2[plm_pos >= lmax] = 0
        plm_pos[plm_pos >= lmax] = lmax-1

        tmp = np.concatenate([np.arange(4,4*nside,4), [4 * nside]*(2*nside-1), np.flip(np.arange(4,4*nside,4))])
        TJI = np.concatenate([[0], np.cumsum(tmp)])

        theta,phi0 = hp.pix2ang(nside, TJI)

        self.nbuff = 1
        if pol == True or niter > 0:
            self.nbuff = 2
        if pol == True and niter > 0:
            self.nbuff = 4

        if all_buff == True:
            self.nbuff = 4
            
        fastSHT.sht_data_alloc(np.array((self.nside, self.lmax, self.nring, self.nsim, self.nbuff, self.pol), dtype=np.int32))
        fastSHT.sht_set_data(ring_tab[:lmax+1],plm_pos, plm_val1, plm_val2, theta, phi0)
        
    def t2alm(self, maps, alms_in=None, method='old', get_alm=True):
        if(alms_in is None and get_alm == True):
            try:
                alms = numba.cuda.pinned_array((self.nsim, self.lmax+1, self.lmax+1), dtype=np.double, strides=None, order='F')
            except:
                alms = np.ones((self.nsim, self.lmax+1, self.lmax+1), dtype=np.double, order='F')
        else:
            alms = alms_in

        if(get_alm == False):
            if(self.niter == 0):
                fastSHT.t2alm(maps)
            elif(method == 'old'):
                fastSHT.t2alm_iter_old(maps, self.niter)
            else:
                fastSHT.t2alm_iter(maps, self.niter)
            return
        else:
            if(self.niter == 0):
                fastSHT.t2alm(maps, alms)
            elif(method == 'old'):
                fastSHT.t2alm_iter_old(maps, self.niter, alms)
            else:
                fastSHT.t2alm_iter(maps, self.niter, alms)

        if(alms_in is None):
            return alms
        else:
            return 0

    def alm2t(self, alms, maps_in=None):
        if(maps_in is None):
            try:
                maps = numba.cuda.pinned_array((self.npix, self.nsim), dtype=np.double, strides=None, order='F')
            except:
                maps = np.ones((self.npix, self.nsim), dtype=np.double, order='F')
        else:
            maps = maps_in
            
        fastSHT.alm2t( alms, maps)
        if(maps_in is None):
            return maps
        else:
            return 0

    def qu2eb(self, Q, U, almEs_in=None, almBs_in=None, method='old'):
        if(almEs_in is None or almBs_in is None):
            try:
                almEs = numba.cuda.pinned_array((self.nsim, self.lmax+1, self.lmax+1), dtype=np.double, strides=None, order='F')
                almBs = numba.cuda.pinned_array((self.nsim, self.lmax+1, self.lmax+1), dtype=np.double, strides=None, order='F')
            except:
                almEs = np.ones((self.nsim, self.lmax+1, self.lmax+1), order='F')
                almBs = np.ones((self.nsim, self.lmax+1, self.lmax+1), order='F')
        else:
            almEs = almEs_in
            almBs = almBs_in
            
        if(self.niter == 0):
            fastSHT.qu2eb(Q, U, almEs, almBs)
        elif(method == 'old'):
            fastSHT.qu2eb_iter_old(Q, U, almEs, almBs, self.niter, 0)

        if(almEs_in is None or almBs_in is None):            
            return (almEs, almBs)
        else:
            return 0

    def qu2e(self, Q, U, almEs_in=None):
        if(almEs_in is None):
            try:
                almEs = numba.cuda.pinned_array((self.nsim, self.lmax+1, self.lmax+1), dtype=np.double, strides=None, order='F')
            except:
                almEs = np.ones((self.nsim, self.lmax+1, self.lmax+1), order='F')
        else:
            almEs = almEs_in

        fastSHT.qu2e(Q, U, almEs)
        
        if(almEs_in is None):            
            return almEs
        else:
            return 0

    def qu2b(self, Q, U, almBs_in=None):
        if(almBs_in is None):
            try:
                almBs = numba.cuda.pinned_array((self.nsim, self.lmax+1, self.lmax+1), dtype=np.double, strides=None, order='F')
            except:
                almBs = np.ones((self.nsim, self.lmax+1, self.lmax+1), order='F')
        else:
            almBs = almBs_in
            
        fastSHT.qu2b(Q, U, almBs)

        if(almBs_in is None):            
            return almBs
        else:
            return 0

        
    def eb2qu(self, almEs, almBs):
        try:
            Q = numba.cuda.pinned_array((self.npix, self.nsim), dtype=np.double, strides=None, order='F')
            U = numba.cuda.pinned_array((self.npix, self.nsim), dtype=np.double, strides=None, order='F')
        except:
            Q = np.ones((self.npix, self.nsim), order='F')
            U = np.ones((self.npix, self.nsim), order='F')


        fastSHT.eb2qu(almEs, almBs, Q, U)
                
        return (Q, U)

    def fix_eb(self, Q, U, mask):
        vid = (np.arange(len(mask))[mask == 1])
        nv = len(vid)
    
        Bmap = np.empty((self.npix, self.nsim), order='F')
        alms = np.empty((self.nsim, self.lmax+1, self.lmax+1), order='F')
        B_template = np.empty((self.npix, self.nsim), order='F')
        flag = 0

        fastSHT.fix_eb(Q, U, mask, Bmap, alms, B_template, vid, self.niter, nv, flag)
        return Bmap
        
    def convert_alm_healpy(self, alms):
        alms_hp = np.empty((2, int((self.lmax+1) * (self.lmax+2) / 2) , self.nsim), order='F')
        fastSHT.convert_alm_healpy1(alms, alms_hp, self.nsim, self.lmax)

        alms_hp = alms_hp[0,:,:] + 1j * alms_hp[1,:,:]

        return alms_hp


    def get_alm(self, alms_in=None):
        if(alms_in is None):
            try:
                alms = numba.cuda.pinned_array((self.nsim, self.lmax+1, self.lmax+1), dtype=np.double, strides=None, order='F')
            except:
                alms = np.ones((self.nsim, self.lmax+1, self.lmax+1), dtype=np.double, order='F')
        else:
            alms = alms_in

        fastSHT.get_alm(alms)
        if(alms_in is None):
            return alms
        else:
            return 0
