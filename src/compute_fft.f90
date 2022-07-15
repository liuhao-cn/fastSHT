!
! Contains: 
! map2fft: compute FFT (forward)
! fft2map: compute map (backward)
!
include 'mkl_dfti.f90'

! Compute map-FFT with real-FFT and PERM-format storage, and swallow the
! transpose job by setting the output strides.
!
! This program only does the forward FFT.
!
  
#ifdef GPU
subroutine map2fft( map, fft, d_fft)
    use MKL_DFTI  
    use sht_data_module
    use cudafor
    implicit none

    type(DFTI_DESCRIPTOR), POINTER  :: FFT_Handle
    real(8), intent(in)             :: map(0:npix-1, 1:nsim)
    real(8), intent(out)            :: fft(1:nsim, 0:npix-1)
    integer(4)                      :: i, j, nfft, status, strides(0:1)

    integer                         :: stat

    real(8), device            :: d_fft(1:nsim, 0:npix-1)
    
    strides(0) = 0
    strides(1) = nsim

    do i=0,nside-1
        nfft   = 4*(i+1)
        j      = nring-i-1
        ! create/recreate descriptor for: 1D-FFT, multi-task, out-of-place,
        ! real-to-real, PERM-format storage, auto-transpose
        status = DftiCreateDescriptor( FFT_Handle, DFTI_DOUBLE,                 DFTI_REAL, 1, nfft )
        status =         DftiSetValue( FFT_Handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_REAL  )    
        status =         DftiSetValue( FFT_Handle, DFTI_PACKED_FORMAT,          DFTI_PERM_FORMAT   )        
        status =         DftiSetValue( FFT_Handle, DFTI_PLACEMENT,              DFTI_NOT_INPLACE   )
        status =         DftiSetValue( FFT_Handle, DFTI_NUMBER_OF_TRANSFORMS,   int(nsim)          )                
        status =         DftiSetValue( FFT_Handle, DFTI_INPUT_DISTANCE,         int(npix)          )
        status =         DftiSetValue( FFT_Handle, DFTI_OUTPUT_DISTANCE,        1                  )
        status =         DftiSetValue( FFT_Handle, DFTI_OUTPUT_STRIDES,         strides            )
        status = DftiCommitDescriptor( FFT_Handle )
        
        status = DftiComputeForward( FFT_Handle, map(h_i1_arr(i):h_i2_arr(i),1), fft(1:nsim,h_i1_arr(i)) ) 
        status = DftiComputeForward( FFT_Handle, map(h_i1_arr(j):h_i2_arr(j),1), fft(1:nsim,h_i1_arr(j)) ) 

        !stat = cudaStreamCreate(cu_streams(i))
        stat = cudaMemcpyAsync(d_fft(1:nsim,h_i1_arr(i)), fft(1:nsim,h_i1_arr(i)), nsim*nfft, cu_streams(i))
        stat = cudaMemcpyAsync(d_fft(1:nsim,h_i1_arr(j)), fft(1:nsim,h_i1_arr(j)), nsim*nfft, cu_streams(i))
        ! free the handle if the next ring will have different size.
        if (i.lt.nside-1) status = DftiFreeDescriptor( FFT_Handle )
     enddo
     do i=nside,3*nside-2
       status = DftiComputeForward( FFT_Handle, map(h_i1_arr(i):h_i2_arr(i),1), fft(1:nsim,h_i1_arr(i)))
       !stat = cudaStreamCreate(cu_streams(i))
       stat = cudaMemcpyAsync(d_fft(1:nsim,h_i1_arr(i)), fft(1:nsim,h_i1_arr(i)), nsim*nfft, cu_streams(i))
    enddo
    status =   DftiFreeDescriptor( FFT_Handle )
#ifdef GPU
    stat = cudaDeviceSynchronize()
#endif
return
end

#else

subroutine map2fft( map, fft )
    use MKL_DFTI  
    use sht_data_module
    implicit none

    type(DFTI_DESCRIPTOR), POINTER  :: FFT_Handle
    real(8), intent(in)             :: map(0:npix-1, 1:nsim)
    real(8), intent(out)            :: fft(1:nsim, 0:npix-1)
    integer(4)                      :: i, j, nfft, status, strides(0:1)
    
    strides(0) = 0
    strides(1) = nsim

    do i=0,nside-1
        nfft   = 4*(i+1)
        j      = nring-i-1
        ! create/recreate descriptor for: 1D-FFT, multi-task, out-of-place,
        ! real-to-real, PERM-format storage, auto-transpose
        status = DftiCreateDescriptor( FFT_Handle, DFTI_DOUBLE,                 DFTI_REAL, 1, nfft )
        status =         DftiSetValue( FFT_Handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_REAL  )    
        status =         DftiSetValue( FFT_Handle, DFTI_PACKED_FORMAT,          DFTI_PERM_FORMAT   )        
        status =         DftiSetValue( FFT_Handle, DFTI_PLACEMENT,              DFTI_NOT_INPLACE   )
        status =         DftiSetValue( FFT_Handle, DFTI_NUMBER_OF_TRANSFORMS,   int(nsim)          )                
        status =         DftiSetValue( FFT_Handle, DFTI_INPUT_DISTANCE,         int(npix)          )
        status =         DftiSetValue( FFT_Handle, DFTI_OUTPUT_DISTANCE,        1                  )
        status =         DftiSetValue( FFT_Handle, DFTI_OUTPUT_STRIDES,         strides            )
        status = DftiCommitDescriptor( FFT_Handle )
        

        status = DftiComputeForward( FFT_Handle, map(i1_arr(i):i2_arr(i),1), fft(1:nsim,i1_arr(i)) ) 
        status = DftiComputeForward( FFT_Handle, map(i1_arr(j):i2_arr(j),1), fft(1:nsim,i1_arr(j)) )         
        ! free the handle if the next ring will have different size.
        if (i.lt.nside-1) status = DftiFreeDescriptor( FFT_Handle )
     enddo
     do i=nside,3*nside-2
       status = DftiComputeForward( FFT_Handle, map(i1_arr(i):i2_arr(i),1), fft(1:nsim,i1_arr(i)))       
    enddo
    status =   DftiFreeDescriptor( FFT_Handle )
return
end

#endif


! Compute backward map-FFT with real-FFT and PERM-format storage, and swallow
! the transpose job by setting the output strides.
!
! This program only does the backward FFT.
! 
SUBROUTINE fft2map( fft, map )
  Use MKL_DFTI
#ifdef GPU
  use sht_data_module, only: nside, npix, nring, nsim, h_i1_arr, h_i2_arr
#else
  use sht_data_module, only: nside, npix, nring, nsim, i1_arr, i2_arr
#endif
    implicit none

    type(DFTI_DESCRIPTOR), POINTER  :: FFT_Handle
    real(8), intent(in)             :: fft(1:nsim, 0:npix-1)
    real(8), intent(out)            :: map(0:npix-1, 1:nsim)
    integer(4)                      :: i, j, nfft, status, strides(0:1)

    strides(0) = 0
    strides(1) = nsim

    do i=0,nside-1
        nfft   = 4*(i+1)
        j      = nring-i-1
        ! create/recreate descriptor for: 1D-FFT, multi-task, out-of-place,
        ! real-to-real, PERM-format storage, auto-transpose         
        status = DftiCreateDescriptor( FFT_Handle, DFTI_DOUBLE,                 DFTI_REAL, 1, nfft )
        status =         DftiSetValue( FFT_Handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_REAL  )    
        status =         DftiSetValue( FFT_Handle, DFTI_PACKED_FORMAT,          DFTI_PERM_FORMAT   )        
        status =         DftiSetValue( FFT_Handle, DFTI_PLACEMENT,              DFTI_NOT_INPLACE   )
        status =         DftiSetValue( FFT_Handle, DFTI_NUMBER_OF_TRANSFORMS,   int(nsim)          )                
        status =         DftiSetValue( FFT_Handle, DFTI_INPUT_DISTANCE,         1                  )
        status =         DftiSetValue( FFT_Handle, DFTI_INPUT_STRIDES,          strides            )
        status =         DftiSetValue( FFT_Handle, DFTI_OUTPUT_DISTANCE,        int(npix)          ) 
        status = DftiCommitDescriptor( FFT_Handle )

#ifdef GPU
        status = DftiComputeBackward( FFT_Handle, fft(1:nsim,h_i1_arr(i)), map(h_i1_arr(i):h_i2_arr(i),1) )  
        status = DftiComputeBackward( FFT_Handle, fft(1:nsim,h_i1_arr(j)), map(h_i1_arr(j):h_i2_arr(j),1) ) 
#else
        status = DftiComputeBackward( FFT_Handle, fft(1:nsim,i1_arr(i)), map(i1_arr(i):i2_arr(i),1) )  
        status = DftiComputeBackward( FFT_Handle, fft(1:nsim,i1_arr(j)), map(i1_arr(j):i2_arr(j),1) ) 
#endif
        ! free the handle if the next ring will have different size.
        if (i.lt.nside-1) status = DftiFreeDescriptor( FFT_Handle )
    enddo

    do i=nside, 3*nside-2
#ifdef GPU
       status = DftiComputeBackward( FFT_Handle, fft(1:nsim,h_i1_arr(i)), map(h_i1_arr(i):h_i2_arr(i),1) )
#else
       status = DftiComputeBackward( FFT_Handle, fft(1:nsim,i1_arr(i)), map(i1_arr(i):i2_arr(i),1) )
#endif
    enddo

    status =   DftiFreeDescriptor( FFT_Handle )

    return
end







