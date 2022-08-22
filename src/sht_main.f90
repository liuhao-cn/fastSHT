!
! contains:
!
! 12 SHT subroutines + their IDL interfaces:
!    t2alm: temperature, forward
!    alm2t: temperature, backward
!
!    QU2EB: polarization, forward, E and B
!    QU2E:  polarization, forward, only E
!    QU2B:  polarization, forward, only B
!
!    EB2QU: polarization, backward, E and B
!    E2QU:  polarization, backward, only E
!    B2QU:  polarization, backward, only B
!
!    t2alm_iter:     temperature, forward, with iteration, new scheme, faster and much more accurate
!    t2alm_iter_old: temperature, forward, with iteration, old scheme, less optimal but same to HEALPix
!    qu2eb_iter:     polarization, forward, with iteration, new scheme, with EB-leakage correction
!    qu2eb_iter_old: polarization, forward, with iteration, old scheme, less optimal but same to HEALPix
!
! 1 subroutine for EB-leakage + its IDL interface:
!    fix_eb
!
! 2 SHT core subroutines
!    mat2alm
!    alm2mat
!
! 5 accessories
!    set_locations
!    apply_mask_mini
!    copy_with_mask_mini
!    lin_resi_mini
!    lin_resi_mini_nomask




SUBROUTINE t2alm( t, alm )
  use sht_data_module
#ifdef GPU
  use cublas_v2
  use cudafor
#endif
  implicit none
  
    ! Main input/output, alm is saved as a real-value square matrix
    real(8), intent(in)  :: t(0:npix-1, 1:nsim) 
    real(8), intent(out) :: alm(1:nsim, 0:lmax, 0:lmax)
    integer(4) :: beginning, end, tot_beginning, tot_end
    real(8) :: rate
    ! scalar variables
    integer(4)           :: m

#ifdef GPU
    integer  :: stat
#endif
    
#ifdef GPU
    call map2fft( t, buff1, d_buff1)
#else
    call map2fft( t, buff1)    
#endif
     do m=0, lmax
        call set_locations(m)

        call get_plm_recursive(0, m)
#ifdef GPU
        call fft2mat_fast(d_buff1, m)
#else
        call fft2mat_fast(buff1, m)
#endif
#ifdef GPU
        call mat2alm(cu_streams(m), d_alm,  plm0, m, 0, 0.d0)
#else
        call mat2alm(alm,  plm0, m, 0, 0.d0)
#endif
     enddo


#ifdef GPU
     alm = d_alm
#endif

    return
end






SUBROUTINE alm2t( alm, t )
    use sht_data_module
#ifdef GPU
    use cublas_v2
    use cudafor
#endif
    implicit none
  
    ! Main input/output, alm is saved as a real-value square matrix
    real(8), intent(out)  :: t(0:npix-1, 1:nsim) 
    real(8), intent(in) :: alm(1:nsim, 0:lmax, 0:lmax)
    integer(4) :: beginning, end, tot_beginning, tot_end
    real(8) :: rate

    ! scalar variables
    integer(4)            :: m, i
#ifdef GPU
    integer  :: stat
    d_alm = alm
#endif

    do m=0, lmax
        call set_locations(m)

        call get_plm_recursive(0, m)

#ifdef GPU
        call alm2mat( cu_streams(m), d_alm, plm0, m, 0, 0.d0)
#else
        call alm2mat( alm, plm0, m, 0, 0.d0)
#endif

#ifdef GPU
        call mat2fft_fast(d_buff1, m)
#else
        call mat2fft_fast(buff1, m)
#endif

    enddo


    ! explicitly set the unused part of buff1 to zero
    !$acc kernels loop
    do i=0,nring-1
#ifdef GPU
       if(2*(lmax+1) .lt. h_nfft_arr(i)) d_buff1(:,2*(lmax+1)+h_i1_arr(i):h_i2_arr(i)) = 0
#else
       if(2*(lmax+1) .lt. nfft_arr(i)) buff1(:,2*(lmax+1)+i1_arr(i):i2_arr(i)) = 0
#endif
    enddo
    !$acc end kernels

#ifdef GPU
    buff1 = d_buff1
#endif
    
    call fft2map( buff1, t )
    return
end














! Similar to t2alm_iter, but use "immediate update" to improve the accuracy.
SUBROUTINE t2alm_iter( t, alm, niter )
    use sht_data_module
#ifdef GPU
    use cublas_v2
    use cudafor
#endif
    implicit none

    ! Main input/output, alm is saved as a real-value square matrix
    real(8), intent(in)     :: t(0:npix-1, 1:nsim)
    real(8), intent(out)    :: alm(1:nsim, 0:lmax, 0:lmax)

    ! scalar variables
    integer(4), intent(in)  :: niter
    integer(4)              :: m, i

#ifdef GPU
    integer  :: stat
#endif

#ifdef GPU
    call map2fft( t, buff1, d_buff1)
#else
    call map2fft( t, buff1 )
#endif

    ! compute initial alm
    do m=0, lmax
        call set_locations(m)
        
        call get_plm_recursive(0, m)
#ifdef GPU

        call fft2mat_fast(d_buff1, m)
#else
        call fft2mat_fast(buff1, m)
#endif
#ifdef GPU

        call mat2alm(cu_streams(m), d_alm,  plm0, m, 0, 0.d0)
#else
        call mat2alm(alm,  plm0, m, 0, 0.d0)
#endif
    enddo

    ! start iteration
    do i=1,niter
        
        do m=0, lmax
            call set_locations(m)
            
            call get_plm_recursive(0, m)

#ifdef GPU
            call alm2mat(cu_streams(m), d_alm, plm0, m, 0, 0.d0)
#else
            call alm2mat( alm, plm0, m, 0, 0.d0)
#endif
#ifdef GPU
            call mat2fft_fast_iter(d_buff1, d_buff2, m)
#else
            call mat2fft_fast_iter(buff1, buff2, m)
#endif
        !-----------------------------------------------------------
        ! The new and old iteration codes are different by the following 4 lines
        ! (new = without them, old = with them)
        !-----------------------------------------------------------
        ! enddo
        ! do m=0, lmax
        !     call set_locations(m)
        !     call get_plm_recursive(0, m)
#ifdef GPU
            call fft2mat_fast(d_buff2, m)            
#else
            call fft2mat_fast(buff2, m)
#endif

#ifdef GPU
            call mat2alm(cu_streams(m), d_alm,  plm0, m, 0, 1.d0)
#else
            call mat2alm(alm,  plm0, m, 0, 1.d0)
#endif
        enddo
    enddo
#ifdef GPU
    alm = d_alm
#endif

    return 
end







! Note: Before HEALPix 3.50, there was a bug that affects the iterations with
! mask, which was fixed in HEALPix 3.50. This programs is consistent to the
! fixed HEALPix routines, so it might give slightly different result to
! HERALPIx 3.31 or earlier ones
SUBROUTINE t2alm_iter_old( t, alm, niter )
  use sht_data_module
#ifdef GPU
  use cublas_v2
  use cudafor
#endif
    implicit none

    ! Main input/output, alm is saved as a real-value square matrix
    real(8), intent(in)     :: t(0:npix-1, 1:nsim)
    real(8), intent(out)    :: alm(1:nsim, 0:lmax, 0:lmax)

    ! scalar variables
    integer(4), intent(in)  :: niter
    integer(4)              :: m, i

#ifdef GPU
    integer  :: stat
#endif

#ifdef GPU
    call map2fft( t, buff1, d_buff1)
#else
    call map2fft( t, buff1 )
#endif
    ! compute initial alm
    do m=0, lmax
        call set_locations(m)

         call get_plm_recursive(0, m)

#ifdef GPU

        call fft2mat_fast(d_buff1, m)
#else
        call fft2mat_fast(buff1, m)
#endif
#ifdef GPU

        call mat2alm(cu_streams(m), d_alm,  plm0, m, 0, 0.d0)
#else
        call mat2alm(alm,  plm0, m, 0, 0.d0)
#endif
    enddo

    ! start iteration
    do i=1,niter
        do m=0, lmax
            call set_locations(m)
            
            call get_plm_recursive(0, m)

#ifdef GPU
            call alm2mat(cu_streams(m), d_alm, plm0, m, 0, 0.d0)
#else
            call alm2mat( alm, plm0, m, 0, 0.d0)
#endif
#ifdef GPU
            call mat2fft_fast_iter(d_buff1, d_buff2, m)
#else
            call mat2fft_fast_iter(buff1, buff2, m)
#endif
         enddo
#ifdef GPU
         stat = cudaDeviceSynchronize()
#endif
        do m=0, lmax
            call set_locations(m)
            
            call get_plm_recursive(0, m)
#ifdef GPU
            call fft2mat_fast(d_buff2, m)            
#else
            call fft2mat_fast(buff2, m)
#endif
            
#ifdef GPU
            call mat2alm(cu_streams(m), d_alm,  plm0, m, 0, 1.d0)
#else
            call mat2alm(alm,  plm0, m, 0, 1.d0)
#endif
        enddo
     enddo
#ifdef GPU
    alm = d_alm
#endif
    return 
end















SUBROUTINE QU2EB( Q, U, almE, almB)
  use sht_data_module
#ifdef GPU
  use cublas_v2
  use cudafor
#endif
    implicit none
    ! Main input/output, alm is saved as real-value square matrix with the new scheme (see fig.1)
    real(8), intent(in)     :: Q(0:npix-1, 1:nsim)
    real(8), intent(in)     :: U(0:npix-1, 1:nsim)
    real(8), intent(out)    :: almE(1:nsim, 0:lmax, 0:lmax)
    real(8), intent(out)    :: almB(1:nsim, 0:lmax, 0:lmax)

    integer(4)              :: m

#ifdef GPU
    integer  :: stat
#endif

#ifdef GPU
    call map2fft(Q, buff1, d_buff1)
    call map2fft(U, buff1, d_buff2)
#else
    call map2fft(Q, buff1)
    call map2fft(U, buff2)
#endif
    do m=0, lmax
        call set_locations(m)
        
        call get_plm_recursive(1, m)
#ifdef GPU
        call fft2mat_fast(d_buff1, m)
        call mat2alm( cu_streams(m), d_alm, plm1, m, 1, 0.d0)
        call mat2alm( cu_streams(m), d_alm2, plm2, m, 3, 0.d0)

        call fft2mat_fast(d_buff2, m)
        call mat2alm( cu_streams(m), d_alm, plm2, m, 2, 1.d0)
        call mat2alm( cu_streams(m), d_alm2, plm1, m, 1, 1.d0)
#else
        call fft2mat_fast(buff1, m)
        call mat2alm( almE, plm1, m, 1, 0.d0)
        call mat2alm( almB, plm2, m, 3, 0.d0)

        call fft2mat_fast(buff2, m)
        call mat2alm( almE, plm2, m, 2, 1.d0)
        call mat2alm( almB, plm1, m, 1, 1.d0)
#endif
    enddo
#ifdef GPU
    almE = d_alm
    almB = d_alm2
#endif
    return 
end
















! !------------------------------------------------
! ! Solve almE and almB iteratively with an EB-leakage correction.
! !
! ! Note: if one skips the EB-leakage correction, then the result of this
! ! program is identical to "qu2eb_iter_old"
! !------------------------------------------------
! SUBROUTINE qu2eb_iter( Q, U, almE, almB, niter )
!     !DIR$ ATTRIBUTES FORCEINLINE :: set_locations
!     use sht_data_module, only: nside, npix, lmax, nsim, nring, plm0, plm1, plm2, &
!         buff1, buff2, buff3, buff4, i1_arr, i2_arr, nfft_arr
!     implicit none

!     ! Main input/output, alm is saved as real-value square matrix with the new scheme (see fig.1)
!     real(8), intent(in)     :: Q(0:npix-1, 1:nsim)
!     real(8), intent(in)     :: U(0:npix-1, 1:nsim)
!     real(8), intent(out)    :: almE(1:nsim, 0:lmax, 0:lmax)
!     real(8), intent(out)    :: almB(1:nsim, 0:lmax, 0:lmax)

!     ! scalar variables
!     integer(4), intent(in)  :: niter

!     real(8), allocatable    :: alm_buff1(:,:,:), alm_buff2(:,:,:), a1_arr(:)
!     integer(4)              :: m, i, inv

!     allocate( alm_buff1(1:nsim, 0:lmax, 0:lmax) ); alm_buff1 = 0
!     allocate( alm_buff2(1:nsim, 0:lmax, 0:lmax) ); alm_buff2 = 0
!     allocate( a1_arr(1:nsim) ); a1_arr = 0

!     ! first use the traditional iteration to compute almE and almB
!     call qu2eb_iter_old(Q, U, almE, almB, niter, 0)

!     ! Then fix the EB-leakage, first E-to-QU. Note: to save time, the results are
!     ! in form of FFT-buffs, not maps

!     do m=0, lmax
!         call set_locations(m)
        
!         call get_plm_recursive(1, m)

!         call alm2mat( almE, plm1, m, 1, 0.d0 )
!         call mat2fft_fast(buff1, m)

!         call alm2mat( almE, plm2, m, 2, 0.d0 )
!         call mat2fft_fast(buff2, m)
!     enddo

!     ! normalize the FFT-buffs so they are suitable as iterative inputs
!     do i=0,nring-1
!         buff1(:,i1_arr(i):i2_arr(i)) = buff1(:,i1_arr(i):i2_arr(i))*nfft_arr(i)
!         buff2(:,i1_arr(i):i2_arr(i)) = buff2(:,i1_arr(i):i2_arr(i))*nfft_arr(i)
!     enddo

!     ! Then QU-to-EB again, with the input being E-mode only. Set skip_fft=1 to
!     ! start directly from the FFT-buff, so Q, U are actually unused. 
!     !
!     ! The result is the harmonic domain template of EB-leakage (alm_buff2)
!     call qu2eb_iter_old(Q, U, alm_buff1, alm_buff2, niter, 1)

!     ! Convert almB (with leakage) and the harmonic domain leakage template
!     ! (alm_buff2) to pixel domain. Use modified form of "alm2t" so they are done
!     ! in one loop.

!     do m=0, lmax
!         call set_locations(m)

!         call get_plm_recursive(0, m)
        
!         call alm2mat( almB,      plm0, m, 0, 0.d0)
!         call mat2fft_fast(buff1, m)

!         call alm2mat( alm_buff2, plm0, m, 0, 0.d0)
!         call mat2fft_fast(buff2, m)
!     enddo

!     call fft2map( buff1, buff3) ! buff3 = B-map with leakage
!     call fft2map( buff2, buff4) ! buff4 = template of leakage

!     ! compute the linear fitting coefficients (a1_arr)
!     call lin_resi_mini_nomask(buff4, buff3, npix, nsim, a1_arr)

!     ! subtract the template from the B-map with coefficients
!     do i=1,nsim
!         almB(i,:,:) = almB(i,:,:) - alm_buff2(i,:,:)*a1_arr(i)
!     enddo

!     deallocate( alm_buff1, alm_buff2, a1_arr )

!     return 
! end

















! ! Note: Before HEALPix 3.50, there was a bug that affects the iterations with
! ! mask, which was fixed in HEALPix 3.50. This programs is consistent to the
! ! fixed HEALPix routines, so it might give slightly different result to
! ! HERALPIx 3.31 or earlier ones
SUBROUTINE qu2eb_iter_old( Q, U, almE, almB, niter, skip_fft )
#ifdef GPU
  use cublas_v2
  use cudafor
#endif

    use sht_data_module
    implicit none

    ! Main input/output, alm is saved as real-value square matrix with the new scheme (see fig.1)
    real(8), intent(in)     :: Q(0:npix-1, 1:nsim)
    real(8), intent(in)     :: U(0:npix-1, 1:nsim)
    real(8), intent(out)    :: almE(1:nsim, 0:lmax, 0:lmax)
    real(8), intent(out)    :: almB(1:nsim, 0:lmax, 0:lmax)

    ! scalar variables
    integer(4), intent(in)  :: niter, skip_fft
    integer(4)              :: m, i, inv
#ifdef GPU
    integer  :: stat
#endif


    if (skip_fft.eq.0) then
#ifdef GPU
       call map2fft(Q, buff1, d_buff1)
       call map2fft(U, buff1, d_buff2)
#else
       call map2fft(Q, buff1)
       call map2fft(U, buff2)
#endif
    endif

    do m=0, lmax
        call set_locations(m)
        
        call get_plm_recursive(1, m)
#ifdef GPU
        call fft2mat_fast(d_buff1, m)
        call mat2alm( cu_streams(m), d_alm, plm1, m, 1, 0.d0)
        call mat2alm( cu_streams(m), d_alm2, plm2, m, 3, 0.d0)

        call fft2mat_fast(d_buff2, m)
        call mat2alm( cu_streams(m), d_alm, plm2, m, 2, 1.d0)
        call mat2alm( cu_streams(m), d_alm2, plm1, m, 1, 1.d0)
#else
        call fft2mat_fast(buff1, m)
        call mat2alm( almE, plm1, m, 1, 0.d0)
        call mat2alm( almB, plm2, m, 3, 0.d0)

        call fft2mat_fast(buff2, m)
        call mat2alm( almE, plm2, m, 2, 1.d0)
        call mat2alm( almB, plm1, m, 1, 1.d0)
#endif
    enddo

    ! start iteration
    do i=1, niter
        do m=0, lmax
            call set_locations(m)
            
            call get_plm_recursive(1, m)

#ifdef GPU
            call alm2mat( cu_streams(m), d_alm, plm1, m, 1, 0.d0)
            call alm2mat( cu_streams(m), d_alm2, plm2, m, 3, 1.d0)
            call mat2fft_fast_iter(d_buff1, d_buff3, m)

            call alm2mat( cu_streams(m), d_alm, plm2, m, 2, 0.d0)
            call alm2mat( cu_streams(m), d_alm2, plm1, m, 1, 1.d0)
            call mat2fft_fast_iter(d_buff2, d_buff4, m)            
#else
            call alm2mat( almE, plm1, m, 1, 0.d0)
            call alm2mat( almB, plm2, m, 3, 1.d0)
            call mat2fft_fast_iter(buff1, buff3, m)

            call alm2mat( almE, plm2, m, 2, 0.d0)
            call alm2mat( almB, plm1, m, 1, 1.d0)
            call mat2fft_fast_iter(buff2, buff4, m)
#endif
        enddo

        do m=0, lmax
            call set_locations(m)
            
            call get_plm_recursive(1, m)
#ifdef GPU
            call fft2mat_fast(d_buff3, m)
            call mat2alm( cu_streams(m), d_alm, plm1, m, 1, 1.d0)
            call mat2alm( cu_streams(m), d_alm2, plm2, m, 3, 1.d0)

            call fft2mat_fast(d_buff4, m)
            call mat2alm( cu_streams(m), d_alm, plm2, m, 2, 1.d0)
            call mat2alm( cu_streams(m), d_alm2, plm1, m, 1, 1.d0)            
#else
            call fft2mat_fast(buff3, m)
            call mat2alm( almE, plm1, m, 1, 1.d0)
            call mat2alm( almB, plm2, m, 3, 1.d0)

            call fft2mat_fast(buff4, m)
            call mat2alm( almE, plm2, m, 2, 1.d0)
            call mat2alm( almB, plm1, m, 1, 1.d0)
#endif
        enddo
    enddo
#ifdef GPU
    almE = d_alm
    almB = d_alm2
#endif
    return 
end

#ifdef GPU
SUBROUTINE device_qu2eb_iter_old( Q, U, almE, almB, niter, skip_fft )
#ifdef GPU
  use cublas_v2
  use cudafor
#endif

    use sht_data_module
    implicit none

    ! Main input/output, alm is saved as real-value square matrix with the new scheme (see fig.1)
    real(8), intent(in)     :: Q(0:npix-1, 1:nsim)
    real(8), intent(in)     :: U(0:npix-1, 1:nsim)
    real(8), intent(out), device    :: almE(1:nsim, 0:lmax, 0:lmax)
    real(8), intent(out), device    :: almB(1:nsim, 0:lmax, 0:lmax)

    ! scalar variables
    integer(4), intent(in)  :: niter, skip_fft
    integer(4)              :: m, i, inv
#ifdef GPU
    integer  :: stat
#endif


    if (skip_fft.eq.0) then
#ifdef GPU
       call map2fft(Q, buff1, d_buff1)
       call map2fft(U, buff1, d_buff2)
#else
       call map2fft(Q, buff1)
       call map2fft(U, buff2)
#endif
    endif

    do m=0, lmax
        call set_locations(m)
        
        call get_plm_recursive(1, m)
#ifdef GPU
        call fft2mat_fast(d_buff1, m)
        call mat2alm( cu_streams(m), almE, plm1, m, 1, 0.d0)
        call mat2alm( cu_streams(m), almB, plm2, m, 3, 0.d0)

        call fft2mat_fast(d_buff2, m)
        call mat2alm( cu_streams(m), almE, plm2, m, 2, 1.d0)
        call mat2alm( cu_streams(m), almB, plm1, m, 1, 1.d0)
#else
        call fft2mat_fast(buff1, m)
        call mat2alm( almE, plm1, m, 1, 0.d0)
        call mat2alm( almB, plm2, m, 3, 0.d0)

        call fft2mat_fast(buff2, m)
        call mat2alm( almE, plm2, m, 2, 1.d0)
        call mat2alm( almB, plm1, m, 1, 1.d0)
#endif
    enddo

    ! start iteration
    do i=1, niter
        do m=0, lmax
            call set_locations(m)
            
            call get_plm_recursive(1, m)

#ifdef GPU
            call alm2mat( cu_streams(m), almE, plm1, m, 1, 0.d0)
            call alm2mat( cu_streams(m), almB, plm2, m, 3, 1.d0)
            call mat2fft_fast_iter(d_buff1, d_buff3, m)

            call alm2mat( cu_streams(m), almE, plm2, m, 2, 0.d0)
            call alm2mat( cu_streams(m), almB, plm1, m, 1, 1.d0)
            call mat2fft_fast_iter(d_buff2, d_buff4, m)            
#else
            call alm2mat( almE, plm1, m, 1, 0.d0)
            call alm2mat( almB, plm2, m, 3, 1.d0)
            call mat2fft_fast_iter(buff1, buff3, m)

            call alm2mat( almE, plm2, m, 2, 0.d0)
            call alm2mat( almB, plm1, m, 1, 1.d0)
            call mat2fft_fast_iter(buff2, buff4, m)
#endif
        enddo

        do m=0, lmax
            call set_locations(m)
            
            call get_plm_recursive(1, m)
#ifdef GPU
            call fft2mat_fast(d_buff3, m)
            call mat2alm( cu_streams(m), almE, plm1, m, 1, 1.d0)
            call mat2alm( cu_streams(m), almB, plm2, m, 3, 1.d0)

            call fft2mat_fast(d_buff4, m)
            call mat2alm( cu_streams(m), almE, plm2, m, 2, 1.d0)
            call mat2alm( cu_streams(m), almB, plm1, m, 1, 1.d0)            
#else
            call fft2mat_fast(buff3, m)
            call mat2alm( almE, plm1, m, 1, 1.d0)
            call mat2alm( almB, plm2, m, 3, 1.d0)

            call fft2mat_fast(buff4, m)
            call mat2alm( almE, plm2, m, 2, 1.d0)
            call mat2alm( almB, plm1, m, 1, 1.d0)
#endif
        enddo
    enddo
    return 
end
#endif











!------------------------------------------------
! fix the EB-leakage when there is a sky mask
! IMPORTANT: to save memory for the GPU code,
! Q, U will be masked after passing to this routine 
!------------------------------------------------

SUBROUTINE fix_eb( Q, U, mask, Bmap, almB, Btpl, vid, niter, nv, flag)
    use sht_data_module
#ifdef GPU
    use cudafor
#endif
    implicit none

    ! Main input/output, alm is saved as real-value square matrix with the new scheme (see fig.1)
    real(8), dimension(0:npix-1, 1:nsim)        :: Q, U, Bmap, Btpl
    real(8), dimension(0:npix-1)                :: mask
    real(8), dimension(1:nsim, 0:lmax, 0:lmax)  :: almB

    ! scalar variables
    integer(4), intent(in)                      :: niter, nv, flag, vid(0:nv-1)
#ifdef GPU
    real(8), allocatable, dimension(:,:,:), device      ::  almE
    integer                                      :: stat
#else
    real(8), allocatable, dimension(:,:,:)      :: almE, alm_buff1
#endif
    real(8), allocatable, dimension(:)          :: a0_arr, a1_arr
    integer(4)                                  :: m, i, inv
    real(8)                                     :: tmp1, tmp2

    ! note: unlike "qu2eb_iter", here we need only one alm-buff, because we don't
    ! need almE in the end, and it is used as buff at certain places.

#ifdef GPU
    allocate( almE     (1:nsim, 0:lmax, 0:lmax) ); almE = 0
    allocate( buff2(1:nsim, 0:npix-1) ); buff2 = 0
    allocate( buff3(1:nsim, 0:npix-1) ); buff3 = 0
    allocate( buff4(1:nsim, 0:npix-1) ); buff4 = 0
#else 
    allocate( almE     (1:nsim, 0:lmax, 0:lmax) ); almE = 0
    allocate( alm_buff1(1:nsim, 0:lmax, 0:lmax) ); alm_buff1 = 0
#endif
    allocate( a0_arr(1:nsim) ); a0_arr = 0
    allocate( a1_arr(1:nsim) ); a1_arr = 0

    Bmap = 0
    almB = 0
    Btpl = 0
    
    ! First copy Q, U to buff3, buff4 with the masked region set to zero. Do
    ! masking on buff3 and buff4 so the input Q, U will not change. Note: although
    ! the operation is very simple, it should be done by a subroutine because
    ! buff3 and buff4 have transposed shapes (same situation in "apply_mask_mini")
    call copy_with_mask_mini(Q, Q, mask, npix, nsim)
    call copy_with_mask_mini(U, U, mask, npix, nsim)

    ! compute almE, almB with the traditional iterations. Note that buff3 and
    ! buff4 serves not only as inputs, but also as computational buffs.
#ifdef GPU
    call device_qu2eb_iter_old(Q, U, almE, d_alm, niter, 0)
    tmp1 = almE(2, 3, 4)
    tmp2 = d_alm(2, 3, 4)
#else
    call qu2eb_iter_old(Q, U, almE, almB, niter, 0)
    tmp1 = almE(2, 3, 4)
    tmp2 = almB(2, 3, 4)
#endif
    

    ! Then E-to-QU, use buff3 and buff4 as the target buffs. Note: after this,
    ! almE can be used as an alm-buff.
    call e2qu( almE, buff3, buff4 )

    ! Apply mask again. Note: unlike "qu2eb_iter", here we must go back to pixel
    ! domain due to the need to apply the mask.
    call apply_mask_mini(buff3, mask, npix, nsim)
    call apply_mask_mini(buff4, mask, npix, nsim)

    ! Then QU-to-EB again with iteration. Again, buff3 and buff4 serves not only
    ! as inputs, but also as computational buffs. The result is the harmonic
    ! domain template of EB-leakage (alm_buff1).
#ifdef GPU
    call device_qu2eb_iter_old(buff3, buff4, almE, d_alm2, niter, 0)
#else
    call qu2eb_iter_old(buff3, buff4, almE, alm_buff1, niter, 0)
#endif
    ! convert almB and the harmonic domain leakage template to pixel domain, use
    ! modified form of "alm2t" so they are done in one loop.

    do m=0, lmax
        call set_locations(m)

        call get_plm_recursive(0, m)

#ifdef GPU
        call alm2mat( cu_streams(m), d_alm, plm0, m, 0, 0.d0)
        call mat2fft_fast(d_buff1, m)

        call alm2mat( cu_streams(m), d_alm2, plm0, m, 0, 0.d0)
        call mat2fft_fast(d_buff2, m)        
#else
        call alm2mat( almB,      plm0, m, 0, 0.d0)
        call mat2fft_fast(buff1, m)

        call alm2mat( alm_buff1, plm0, m, 0, 0.d0)
        call mat2fft_fast(buff2, m)
#endif
    enddo

    !$acc kernels loop
    do i=0,nring-1
#ifdef GPU
        if(2*(lmax+1) .lt. h_nfft_arr(i)) d_buff1(:,2*(lmax+1)+h_i1_arr(i):h_i2_arr(i)) = 0
        if(2*(lmax+1) .lt. h_nfft_arr(i)) d_buff2(:,2*(lmax+1)+h_i1_arr(i):h_i2_arr(i)) = 0
#else
        if(2*(lmax+1) .lt. nfft_arr(i)) buff1(:,2*(lmax+1)+i1_arr(i):i2_arr(i)) = 0
        if(2*(lmax+1) .lt. nfft_arr(i)) buff2(:,2*(lmax+1)+i1_arr(i):i2_arr(i)) = 0
#endif
     enddo
     !$end kernels

    ! if(lmax.lt.2*nside) forall(i=0:npix-1, acc_flag1(i).eq.0) buff1(:,i) = 0
    ! if(lmax.lt.2*nside) forall(i=0:npix-1, acc_flag2(i).eq.0) buff2(:,i) = 0

#ifdef GPU
     buff1 = d_buff1
     buff2 = d_buff2
#endif
     
    call fft2map( buff1, buff3) ! buff3 = B-map with leakage
    call fft2map( buff2, Btpl)  ! buff4 = template of leakage

    ! final linear fitting, only for the available region
    call lin_resi_mini(Btpl, buff3, Bmap, npix, nsim, vid, nv, a0_arr, a1_arr)

    ! compute the final almB only when flag=1
    ! if (flag.eq.1) then
    !     do i=1,nsim
    !         almB(i,:,:) = almB(i,:,:) - alm_buff1(i,:,:)*a1_arr(i)
    !         almB(i,0,0) = almB(i,0,0) - a0_arr(i)*sqrt(4.d0*pi)
    !     enddo
    ! endif

    ! deallocate( almE, alm_buff1, a0_arr, a1_arr )

    return 
end













! SUBROUTINE QU2E( Q, U, almE)
!     !DIR$ ATTRIBUTES FORCEINLINE :: set_locations
!     use sht_data_module, only: nside, npix, lmax, nsim, nring, plm1, plm2, buff1, buff2
!     implicit none
!     ! Main input/output, alm is saved as real-value square matrix with the new scheme (see fig.1)
!     real(8), intent(in)     :: Q(0:npix-1, 1:nsim)
!     real(8), intent(in)     :: U(0:npix-1, 1:nsim)
!     real(8), intent(out)    :: almE(1:nsim, 0:lmax, 0:lmax)

!     integer(4) :: m

!     call map2fft(Q, buff1)
!     call map2fft(U, buff2)

!     do m=0, lmax
!         call set_locations(m)
        
!         call get_plm_recursive(1, m)

!         call fft2mat_fast(buff1, m)
!         call mat2alm( almE, plm1, m, 1, 0.d0)

!         call fft2mat_fast(buff2, m)
!         call mat2alm( almE, plm2, m, 2, 1.d0)
!     enddo

!     return 
! end










! SUBROUTINE QU2B( Q, U, almB)
!     !DIR$ ATTRIBUTES FORCEINLINE :: set_locations
!     use sht_data_module, only: nside, npix, lmax, nsim, nring, plm1, plm2, buff1, buff2
!     implicit none
!     ! Main input/output, alm is saved as real-value square matrix with the new scheme (see fig.1)
!     real(8), intent(in)     :: Q(0:npix-1, 1:nsim)
!     real(8), intent(in)     :: U(0:npix-1, 1:nsim)
!     real(8), intent(out)    :: almB(1:nsim, 0:lmax, 0:lmax)

!     integer(4) :: m

!     call map2fft(Q, buff1)
!     call map2fft(U, buff2)

!     do m=0, lmax
!         call set_locations(m)

!         call get_plm_recursive(1, m)

!         call fft2mat_fast(buff2, m)
!         call mat2alm( almB, plm1, m, 1, 0.d0)

!         call fft2mat_fast(buff1, m)
!         call mat2alm( almB, plm2, m, 3, 1.d0)
!     enddo

!     return 
! end














SUBROUTINE EB2QU( almE, almB, Q, U )
    use sht_data_module
#ifdef GPU
    use cublas_v2
    use cudafor
#endif
    implicit none

    ! Main input/output, alm is saved as real-value square matrix with the new scheme (see fig.1)
    real(8), intent(in)    :: almE(1:nsim, 0:lmax, 0:lmax)
    real(8), intent(in)    :: almB(1:nsim, 0:lmax, 0:lmax)
    real(8), intent(out)   :: Q(0:npix-1, 1:nsim)
    real(8), intent(out)   :: U(0:npix-1, 1:nsim)

    integer(4) :: m, i

#ifdef GPU
    integer  :: stat
    d_alm = almE
    d_alm2 = almB
#endif


    do m=0, lmax
        call set_locations(m)
        
        call get_plm_recursive(1, m)

#ifdef GPU
        call alm2mat( cu_streams(m), d_alm, plm1, m, 1, 0.d0)
        call alm2mat( cu_streams(m), d_alm2, plm2, m, 3, 1.d0)
        call mat2fft_fast(d_buff1, m)

        call alm2mat( cu_streams(m), d_alm, plm2, m, 2, 0.d0)
        call alm2mat( cu_streams(m), d_alm2, plm1, m, 1, 1.d0)
        call mat2fft_fast(d_buff2, m)
#else
        call alm2mat( almE, plm1, m, 1, 0.d0)
        call alm2mat( almB, plm2, m, 3, 1.d0)
        call mat2fft_fast(buff1, m)

        call alm2mat( almE, plm2, m, 2, 0.d0)
        call alm2mat( almB, plm1, m, 1, 1.d0)
        call mat2fft_fast(buff2, m)
#endif
    enddo

    !$acc kernels loop
    do i=0,nring-1
#ifdef GPU
        if(2*(lmax+1) .lt. h_nfft_arr(i)) d_buff1(:,2*(lmax+1)+h_i1_arr(i):h_i2_arr(i)) = 0
        if(2*(lmax+1) .lt. h_nfft_arr(i)) d_buff2(:,2*(lmax+1)+h_i1_arr(i):h_i2_arr(i)) = 0       
#else
        if(2*(lmax+1) .lt. nfft_arr(i)) buff1(:,2*(lmax+1)+i1_arr(i):i2_arr(i)) = 0
        if(2*(lmax+1) .lt. nfft_arr(i)) buff2(:,2*(lmax+1)+i1_arr(i):i2_arr(i)) = 0
#endif
     enddo
     !$acc end kernels

#ifdef GPU
     buff1 = d_buff1
     buff2 = d_buff2
#endif
     
     call fft2map( buff1, Q )
     call fft2map( buff2, U )

     return 
end

















SUBROUTINE E2QU( almE, Q, U )
  use sht_data_module
#ifdef GPU
    use cudafor
#endif
    implicit none
    ! Main input/output, alm is saved as real-value square matrix with the new scheme (see fig.1)
#ifdef GPU
    real(8), intent(in), device   :: almE(1:nsim, 0:lmax, 0:lmax)
    real(8), intent(out)   :: Q(0:npix-1, 1:nsim)
    real(8), intent(out)   :: U(0:npix-1, 1:nsim)
    integer  :: stat
#else
    real(8), intent(in)    :: almE(1:nsim, 0:lmax, 0:lmax)
    real(8), intent(out)   :: Q(0:npix-1, 1:nsim)
    real(8), intent(out)   :: U(0:npix-1, 1:nsim)
#endif
    integer(4) :: m, i

    do m=0, lmax

        call set_locations(m)
        
        call get_plm_recursive(1, m)
#ifdef GPU
        call alm2mat( cu_streams(m), almE, plm1, m, 1, 0.d0 )
        call mat2fft_fast(d_buff1, m)

        call alm2mat( cu_streams(m), almE, plm2, m, 2, 0.d0 )
        call mat2fft_fast(d_buff2, m)
#else
        call alm2mat( almE, plm1, m, 1, 0.d0 )
        call mat2fft_fast(buff1, m)

        call alm2mat( almE, plm2, m, 2, 0.d0 )
        call mat2fft_fast(buff2, m)
#endif
    enddo

    !$acc kernels loop
    do i=0,nring-1
#ifdef GPU
        if(2*(lmax+1) .lt. h_nfft_arr(i)) d_buff1(:,2*(lmax+1)+h_i1_arr(i):h_i2_arr(i)) = 0
        if(2*(lmax+1) .lt. h_nfft_arr(i)) d_buff2(:,2*(lmax+1)+h_i1_arr(i):h_i2_arr(i)) = 0
#else
        if(2*(lmax+1) .lt. nfft_arr(i)) buff1(:,2*(lmax+1)+i1_arr(i):i2_arr(i)) = 0
        if(2*(lmax+1) .lt. nfft_arr(i)) buff2(:,2*(lmax+1)+i1_arr(i):i2_arr(i)) = 0

#endif
    enddo
    !$end kernels
#ifdef GPU
    buff1 = d_buff1
    buff2 = d_buff2
#endif

    call fft2map( buff1, Q )
    call fft2map( buff2, U )

    return 
end












! SUBROUTINE B2QU( almB, Q, U )
!     !DIR$ ATTRIBUTES FORCEINLINE :: set_locations
!     use sht_data_module, only: nside, npix, lmax, nsim, nring, plm1, plm2, &
!         buff1, buff2, i1_arr, i2_arr, nfft_arr
!     implicit none
!     ! Main input/output, alm is saved as real-value square matrix with the new scheme (see fig.1)
!     real(8), intent(in)    :: almB(1:nsim, 0:lmax, 0:lmax)
!     real(8), intent(out)   :: Q(0:npix-1, 1:nsim)
!     real(8), intent(out)   :: U(0:npix-1, 1:nsim)

!     integer(4) :: m, i

!     do m=0, lmax
!         call set_locations(m)
        
!         call get_plm_recursive(1, m)

!         call alm2mat( almB, plm1, m, 1, 0.d0)
!         call mat2fft_fast(buff2, m)

!         call alm2mat( almB, plm2, m, 3, 0.d0)
!         call mat2fft_fast(buff1, m)
!     enddo

!     do i=0,nring-1
!         if(2*(lmax+1) .lt. nfft_arr(i)) buff1(:,2*(lmax+1)+i1_arr(i):i2_arr(i)) = 0
!         if(2*(lmax+1) .lt. nfft_arr(i)) buff2(:,2*(lmax+1)+i1_arr(i):i2_arr(i)) = 0
!     enddo

!     call fft2map( buff1, Q )
!     call fft2map( buff2, U )

!     return 
! end























