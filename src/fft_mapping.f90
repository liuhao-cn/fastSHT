!
! contains the24 codes for mapping the mapFFT:
!
! fft2mat_fast
! mat2fft_fast
! mat2fft_fast_iter
!
!
!
!------------------------------------------------------------------------- 
! parse the Perm format real map-FFT to convert the map-FFT to FFT-matrix
!
! use precomputed table and factor to work faster
!
! use the precomputed tables and factors to work faster
!-------------------------------------------------------------------------


subroutine fft2mat_fast(mapfft, m)
  use sht_data_module, only: nring, nsim, npix, i0, mat_cache1, mat_cache2, tab, fac1, fac2, fac3, fac4
#ifdef GPU
  use cudafor
#endif
    implicit none
#ifdef GPU
    real(8),    dimension(1:nsim, 0:npix-1), device, intent(in)  :: mapfft
#else
    real(8),    dimension(1:nsim, 0:npix-1),  intent(in)  :: mapfft
#endif

    integer(4), intent(in)  :: m
    integer(4)              :: i, j, k

#ifdef GPU
    i = cudaDeviceSynchronize()
#endif

    !$acc kernels
    !$acc loop private(j)
#ifndef GPU
    !$omp parallel do private(j)
#endif
    do k=1,nsim
    do i=i0,(nring+1)/2-1
        j = tab(i,m)
            mat_cache1(k,i) = mapfft(k,j)*fac1(i,m) + mapfft(k,j+1)*fac2(i,m)
            mat_cache2(k,i) = mapfft(k,j)*fac3(i,m) + mapfft(k,j+1)*fac4(i,m)
        enddo
     enddo
#ifndef GPU
    !$omp end parallel do
#endif
     !$acc loop private(j)
#ifndef GPU
    !$omp parallel do private(j)
#endif
     do k=1,nsim
        do i=(nring+1)/2,nring-1-i0
           j = tab(i,m)
           mat_cache1(k,i) = mapfft(k,j)*fac1(nring-1-i,m) + mapfft(k,j+1)*fac2(nring-1-i,m)
           mat_cache2(k,i) = mapfft(k,j)*fac3(nring-1-i,m) + mapfft(k,j+1)*fac4(nring-1-i,m)
        enddo
     enddo
#ifndef GPU
    !$omp end parallel do
#endif
    !$acc end kernels
end subroutine








!-------------------------------------------------------------------------
! parse the Perm format real map-FFT to save the FFT-matrix into the map-FFT
! use precomputed table and factor to work faster 
!
! for the backward transform without iteration
! written in a way that the zeroing is need only on a small array
!-------------------------------------------------------------------------
subroutine mat2fft_fast(mapfft, m)
    use sht_data_module, only: nside, nring, lmax, nsim, npix, i0, mat_cache1, mat_cache2, &
        nfft_arr, tabI, fac1I, fac2I, fac3I, fac4I
#ifdef GPU
    use cudafor
#endif
    implicit none

#ifdef GPU
    real(8), device,   intent(out) :: mapfft(1:nsim, 0:npix-1)
#else
    real(8),    intent(out) :: mapfft(1:nsim, 0:npix-1)
#endif
    integer(4), intent(in)  :: m
    integer(4)              :: i, j, jf, k, stat
    real(8)                 :: fac

#ifdef GPU
    stat = cudaDeviceSynchronize()
#endif

    !$acc kernels
    !$acc loop 
#ifndef GPU
    !$omp parallel do private(j,jf,fac)
#endif
    do k=1,nsim
       !$acc loop  independent private(j,jf,fac)
       do i=i0,nring-1-i0
          j = tabI(i,m)
          jf = min(i, nring-1-i)
          if(2*(m+1).gt.nfft_arr(i)) then
             fac = 1
          else
             fac = 0
          endif
          mapfft(k,  j) = mapfft(k,  j)*fac + mat_cache1(k,i)*fac1I(jf,m) + mat_cache2(k,i)*fac2I(jf,m)
          mapfft(k,1+j) = mapfft(k,1+j)*fac + mat_cache1(k,i)*fac3I(jf,m) + mat_cache2(k,i)*fac4I(jf,m)
       enddo
    enddo
#ifndef GPU
    !$omp end parallel do
#endif
    !$acc end kernels
end subroutine





!-------------------------------------------------------------------------
! parse the Perm format real map-FFT to save the FFT-matrix into the map-FFT
! use precomputed table and factor to work faster 
!
! for the backward transform in iterations
! written in a way that the zeroing is need only on a small array
!-------------------------------------------------------------------------

subroutine mat2fft_fast_iter(mapfft1,mapfft2, m)
    use sht_data_module, only: nside, nring, lmax, nsim, npix, i0, &
         mat_cache1, mat_cache2, nfft_arr, tabI, fac1I, fac2I, fac3I, fac4I
#ifdef GPU
    use cudafor
#endif
    implicit none

#ifdef GPU
     real(8),    dimension(1:nsim, 0:npix-1), device  :: mapfft1
     real(8),    dimension(1:nsim, 0:npix-1), device :: mapfft2
#else
     real(8),    dimension(1:nsim, 0:npix-1), intent(in)  :: mapfft1
     real(8),    dimension(1:nsim, 0:npix-1), intent(out) :: mapfft2
#endif
    
    integer(4), intent(in)  :: m
    integer(4)              :: i, k, stat
    real(8)                 :: fac, tmp
    integer(4)             :: j, jf, nfft
#ifdef GPU
    stat = cudaDeviceSynchronize()
#endif

    !$acc kernels
    !$acc loop 
#ifndef GPU
    !$omp parallel do private(j,jf,nfft)
#endif
    do k=1,nsim
       !$acc loop private(j,jf,nfft) independent
       do i=i0,nring-1-i0
          j = tabI(i,m)
          jf = min(i, nring-1-i)
          nfft = nfft_arr(i)
          if(2*(m+1).gt.nfft) then
             mapfft2(k,  j) = mapfft2(k, j) - mat_cache1(k,i)*fac1I(jf,m)*nfft - mat_cache2(k,i)*fac2I(jf,m)*nfft
             mapfft2(k,1+j) = mapfft2(k,1+j) - mat_cache1(k,i)*fac3I(jf,m)*nfft - mat_cache2(k,i)*fac4I(jf,m)*nfft
          else
             mapfft2(k,  j) = mapfft1(k,  j) - mat_cache1(k,i)*fac1I(jf,m)*nfft - mat_cache2(k,i)*fac2I(jf,m)*nfft
             mapfft2(k,1+j) = mapfft1(k,1+j) - mat_cache1(k,i)*fac3I(jf,m)*nfft - mat_cache2(k,i)*fac4I(jf,m)*nfft
          endif
       enddo
    enddo
#ifndef GPU
    !$omp end parallel do
#endif
    !$acc end kernels
end subroutine




































































! obsolete codes below


! SUBROUTINE mat2fft(argc, argv)            ! interface to IDL  
! !DEC$ ATTRIBUTES DLLEXPORT::mat2fft
! ! 
! INTEGER*8   argc, argv(*)                 ! Argc and Argv are L64 integers    

! CALL mat2fft_run( %VAL(argv( 1)), %VAL(argv( 2)), %VAL(argv( 3)), &
!                   %VAL(argv( 4)), %VAL(argv( 5)), %VAL(argv( 6)), &
!                   %VAL(argv( 7)), %VAL(argv( 8)), %VAL(argv( 9)) )

! RETURN     
! END    

! ! parse the Perm format real map-FFT to save the FFT-matrix into the map-FFT for m in [m1, m2]. 
! ! flag=0: normal mode. flag=1: iterative mode
! ! 
! ! Now this is obsolete
! !
! subroutine mat2fft_run(mapfft, mat_re, mat_im, npix, nsim, nring, m1, m2, phi0, flag)
! !dir$ attributes forceinline :: get_ring_info
! !DEC$ ATTRIBUTES DLLEXPORT::mat2fft_run
! implicit none

! real(8)     :: mapfft(1:nsim, 0:npix-1)
! real(8)     :: mat_re(1:nsim, 0:nring-1, m1:m2)
! real(8)     :: mat_im(1:nsim, 0:nring-1, m1:m2)
! real(8)     :: cc, ss, phi0(0:nring-1)
! integer(4)  :: npix, nsim, nring, m, m1, m2, m3, pos, i, i1, i2, nside, nfft, fac_m, flag, fac

! nside = nint( (npix/12)**0.5, 4 ) 

! do i=0,nring-1       
!     ! nfft is the actual FFT-size
!     call get_ring_info(nside, i, i1, i2, nfft)
!     if (flag.eq.0) then
!         fac = 1         ! for normal computation
!     else
!         fac = -nfft     ! for iteration, must be used with proper setting of the initial value in "mapfft"
!     endif
!     ! to parse the result, first map m to the unit circle, and then map to the
!     ! upper unit circle because the inputs are real (called conjugate-even)
!     do m=m1,m2
!         ! In most cases, due to the Perm scheme definition, the factor of
!         ! storage is 1. The only exception is when none-zero m is mapped to
!         ! the real-axis. In this case, both m and -m will contribute to the
!         ! result, so the factor of storage is 2, which is not covered by the
!         ! Perm scheme, so we must manually set fac_m=2 in this case.
!         if (m.eq.0) then
!             fac_m = 1
!         else
!             fac_m = 2
!         endif
!         ! "pos" is the mapped position on the unit circle, with 0-->0, nfft/2-->Pi, nfft-->2Pi.
!         pos  = mod( m, nfft)
!         ! constant phase angle attached to the starting point of the ring
!         cc   = cos( m*phi0(i) )
!         ss   = sin( m*phi0(i) )
!         ! the correct position on the upper unit circle (matching the Perm format) is "m3"
!         ! note that "pos" already satisfy pos<nfft.
!         if (pos.eq.0) then                      
!             ! zero frequency       
!             m3 = i1
!             mapfft(:,  m3) = mapfft(:,  m3) + (mat_re(:,i,m)*cc-mat_im(:,i,m)*ss)*fac_m*fac
!         else if (pos.lt.nfft/2) then            
!             ! positive frequency (already in the upper unit circle), +1 for imaginary part
!             m3 = i1+pos*2
!             mapfft(:,  m3) = mapfft(:,  m3) + (mat_re(:,i,m)*cc-mat_im(:,i,m)*ss)*fac
!             mapfft(:,1+m3) = mapfft(:,1+m3) + (mat_re(:,i,m)*ss+mat_im(:,i,m)*cc)*fac 
!         else if (pos.eq.nfft/2) then            
!             ! Nyquist frequency, in the Perm format (even size), this is next to zero frequency
!             m3 = i1+1
!             mapfft(:,  m3) = mapfft(:,  m3) + (mat_re(:,i,m)*cc-mat_im(:,i,m)*ss)*fac_m*fac
!         else                                    
!             ! negative frequency, take conjugate and flip to the upper half of the unit circle
!             m3 = i1+(nfft-pos)*2
!             mapfft(:,  m3) = mapfft(:,  m3) + (mat_re(:,i,m)*cc-mat_im(:,i,m)*ss)*fac
!             mapfft(:,1+m3) = mapfft(:,1+m3) - (mat_re(:,i,m)*ss+mat_im(:,i,m)*cc)*fac
!         endif             
!     enddo
! enddo

! end subroutine




















! !
! ! obsolete code
! !
! ! parse the Perm format real-value map-FFT, load values into a rectangular 
! ! FFT-matrix for m in [m1, m2]. Note that map-FFT is already transposed
! !
! ! NOW THIS IS UNUSED ANYMORE, and is replaced by the routine that used
! ! precomputed tables and factors called "fft2mat_fast" (below)
! !
! subroutine fft2mat_run(mapfft, mat_re, mat_im, npix, nsim, nring, m1, m2, phi0)
! !dir$ attributes forceinline :: get_ring_info
! !DEC$ ATTRIBUTES DLLEXPORT::fft2mat_run
! implicit none

! real(8)     :: mapfft(1:nsim, 0:npix-1)
! real(8)     :: mat_re(1:nsim, 0:nring-1, m1:m2)
! real(8)     :: mat_im(1:nsim, 0:nring-1, m1:m2)
! real(8)     :: cc, ss, phi0(0:nring-1)
! integer(4)  :: npix, nsim, nring, m, m1, m2, pos, i, i1, i2, nside, nfft

! nside = nint( (npix/12)**0.5, 4 ) 

! do i=0,nring-1       
!     call get_ring_info(nside, i, i1, i2, nfft)! nfft is the actual FFT-size
!     do m=m1,m2    
!         pos  = mod(m, nfft)
!         cc   = cos( -m*phi0(i) )
!         ss   = sin( -m*phi0(i) )
!         if (pos.eq.0) then                    ! zero FFT-frequency       
!             mat_re(:,i,m) = mapfft(:,i1)*cc
!             mat_im(:,i,m) = mapfft(:,i1)*ss
!         else if (pos.lt.nfft/2) then          ! positive FFT-frequency
!             mat_re(:,i,m) = mapfft(:,i1+pos*2)*cc - mapfft(:,i1+1+pos*2)*ss
!             mat_im(:,i,m) = mapfft(:,i1+pos*2)*ss + mapfft(:,i1+1+pos*2)*cc
!         else if (pos.eq.nfft/2) then          ! Nyquist FFT-frequency
!             mat_re(:,i,m) = mapfft(:,i1+1)*cc
!             mat_im(:,i,m) = mapfft(:,i1+1)*ss
!         else                                  ! negative FFT-frequency            
!             mat_re(:,i,m) = mapfft(:,i1+(nfft-pos)*2)*cc + mapfft(:,i1+1+(nfft-pos)*2)*ss
!             mat_im(:,i,m) = mapfft(:,i1+(nfft-pos)*2)*ss - mapfft(:,i1+1+(nfft-pos)*2)*cc
!         endif
!     enddo
! enddo

! end subroutine












! subroutine mat2fft_fast(mapfft, acc_flag, m)
!     use sht_data_module, only: nside, nring, lmax, nsim, npix, i0, mat_cache1, mat_cache2, tabI, fac1I, fac2I, fac3I, fac4I
!     !DEC$ ATTRIBUTES DLLEXPORT::mat2fft_fast
!     implicit none

!     real(8),    dimension(1:nsim, 0:npix-1), intent(out) :: mapfft
!     integer(4), dimension(0:npix-1)                      :: acc_flag

!     integer(4), intent(in)  :: m
!     integer(4)              :: i, j

!     do i=i0,(nring+1)/2-1
!         j = tabI(i,m)
!         mapfft(:,j) = mapfft(:,j)*acc_flag(j) + mat_cache1(:,i)*fac1I(i,m) + mat_cache2(:,i)*fac2I(i,m)
!         acc_flag(j) = 1
!         ! if (acc_flag(j).eq.0) then
!         !     ! set the value        
!         !     mapfft(:,j) = mat_cache1(:,i)*fac1I(i,m) + mat_cache2(:,i)*fac2I(i,m)
!         !     acc_flag(j) = 1
!         ! else
!         !     ! accumulate the value
!         !     mapfft(:,j) = mapfft(:,j) + mat_cache1(:,i)*fac1I(i,m) + mat_cache2(:,i)*fac2I(i,m)
!         ! endif

!         mapfft(:,1+j) = mapfft(:,1+j)*acc_flag(1+j) + mat_cache1(:,i)*fac3I(i,m) + mat_cache2(:,i)*fac4I(i,m)
!         acc_flag(1+j) = 1
!         ! if (acc_flag(1+j).eq.0) then      
!         !     mapfft(:,1+j) = mat_cache1(:,i)*fac3I(i,m) + mat_cache2(:,i)*fac4I(i,m)        
!         !     acc_flag(1+j) = 1
!         ! else
!         !     mapfft(:,1+j) = mapfft(:,1+j) + mat_cache1(:,i)*fac3I(i,m) + mat_cache2(:,i)*fac4I(i,m)
!         ! endif
!     enddo

!     do i=(nring+1)/2,nring-1-i0
!         j = tabI(i,m)
!         mapfft(:,j) = mapfft(:,j)*acc_flag(j) + mat_cache1(:,i)*fac1I(nring-1-i,m) + mat_cache2(:,i)*fac2I(nring-1-i,m)
!         acc_flag(j) = 1
!         ! if (acc_flag(j).eq.0) then
!         !     ! set the value        
!         !     mapfft(:,j) = mat_cache1(:,i)*fac1I(nring-1-i,m) + mat_cache2(:,i)*fac2I(nring-1-i,m)
!         !     acc_flag(j) = 1
!         ! else
!         !     ! accumulate the value
!         !     mapfft(:,j) = mapfft(:,j) + mat_cache1(:,i)*fac1I(nring-1-i,m) + mat_cache2(:,i)*fac2I(nring-1-i,m)
!         ! endif
        
!         mapfft(:,1+j) = mapfft(:,1+j)*acc_flag(1+j) + mat_cache1(:,i)*fac3I(nring-1-i,m) + mat_cache2(:,i)*fac4I(nring-1-i,m)            
!         acc_flag(1+j) = 1
!         ! if (acc_flag(1+j).eq.0) then   
!         !     ! set the value       
!         !     mapfft(:,1+j) = mat_cache1(:,i)*fac3I(nring-1-i,m) + mat_cache2(:,i)*fac4I(nring-1-i,m)
!         !     acc_flag(1+j) = 1
!         ! else
!         !     ! accumulate the value
!         !     mapfft(:,1+j) = mapfft(:,1+j) + mat_cache1(:,i)*fac3I(nring-1-i,m) + mat_cache2(:,i)*fac4I(nring-1-i,m)            
!         ! endif
!     enddo

! end subroutine
