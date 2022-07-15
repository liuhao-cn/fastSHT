! 1 SHT core subroutines
!    mat2alm
!    alm2mat
!
! 2 accessories
!    set_locations
!    apply_mask_mini
!    copy_with_mask_mini
!    lin_resi_mini
!    lin_resi_mini_nomask


! ---------------------------------------
! This is the core subprogram of SHT, and it is normally unnecessary to call
! this program from external wrappers.
!
! ---------------------------------------
!
! ---------------------------------------
! convert mapFFT to alm, for only one m.
!
! type: used to control the value of a1, a2 and swap of real/imaginary parts,
! designed for polarization
!
! beta: parameter of dgemm, controls "set" or "accumulate"
!
! note that the content of plm should also be for one m
! ---------------------------------------
#ifdef GPU
subroutine mat2alm(stream, alm, plm, m, type, beta)
    use sht_data_module
    use cudafor
    use cublas_v2
    implicit none

    !type(cublasHandle) :: handle
    integer(cuda_stream_kind):: stream
    ! Main input/output, note that alm is flattened from 3D to 2D.
    real(8), dimension(0:nring-1,-2:lmax),device   :: plm
    real(8), dimension(1:nsim, 0:lmax*(lmax+2)),device :: alm

    ! scalar variables
    integer(4), intent(in)  :: m, type
    real(8),    intent(in)  :: beta
    integer                 :: stat, i, j
    real(8)  :: a1p, a2p

    stat = cublasSetStream(handle, stream)

    if (h_swap(0, type).eq.0) then
        stat = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, nsim, n, rlen, h_a1(0,type), mat_cache1(:,i0), nsim, plm(i0:i1,m), nring, beta, alm(:,p1), nsim)
        if (m.gt.0) &
        stat = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, nsim, n, rlen, h_a2(0,type), mat_cache2(:,i0), nsim, plm(i0:i1,m), nring, beta, alm(:,p2), nsim)
     else
        stat = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, nsim, n, rlen, h_a1(0,type), mat_cache2(:,i0), nsim, plm(i0:i1,m), nring, beta, alm(:,p1), nsim)
        if (m.gt.0) &
        stat = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, nsim, n, rlen, h_a2(0,type), mat_cache1(:,i0), nsim, plm(i0:i1,m), nring, beta, alm(:,p2), nsim)
     endif

end subroutine
#else
SUBROUTINE mat2alm( alm, plm, m, type, beta)
    use sht_data_module, only: nside, nsim, npix, nring, lmax, mat_cache1, mat_cache2, &
                               i0, i1, p1, p2, n, rlen, swap, a1, a2
    implicit none

    ! Main input/output, note that alm is flattened from 3D to 2D.
    real(8), dimension(0:nring-1,-2:lmax)       :: plm
    real(8), dimension(1:nsim, 0:lmax*(lmax+2)) :: alm

    ! scalar variables
    integer(4), intent(in)  :: m, type
    real(8),    intent(in)  :: beta

    ! compute alm using the BLAS function dgemm provided by MKL
    !
    ! IMPORTATN: on exit, the output alm already has the result at the correct
    ! positions. No need to copy anything.
    if (swap(0, type).eq.0) then
        call dgemm("N", "N", nsim, n, rlen, a1(0,type), mat_cache1(:,i0), nsim, plm(i0:i1,m), nring, beta, alm(:,p1), nsim)
        if (m.gt.0) &
        call dgemm("N", "N", nsim, n, rlen, a2(0,type), mat_cache2(:,i0), nsim, plm(i0:i1,m), nring, beta, alm(:,p2), nsim)
    else
        call dgemm("N", "N", nsim, n, rlen, a1(0,type), mat_cache2(:,i0), nsim, plm(i0:i1,m), nring, beta, alm(:,p1), nsim)
        if (m.gt.0) &
        call dgemm("N", "N", nsim, n, rlen, a2(0,type), mat_cache1(:,i0), nsim, plm(i0:i1,m), nring, beta, alm(:,p2), nsim)
    endif    

end subroutine
#endif



! ---------------------------------------
! convert alm to mapFFT, for only one m.
!
! type: used to control the value of a1, a2 and swap of real/imaginary parts.
! 
! beta: parameter of dgemm, controls "set" or "accumulate"
!
! note that the content of plm should also be for one m
! ---------------------------------------
#ifdef GPU
SUBROUTINE alm2mat(stream, alm, plm, m, type, beta)
    use sht_data_module
    use cudafor
    use cublas_v2
    implicit none

    integer(cuda_stream_kind):: stream
    ! Main input/output, note that alm is flattened from 3D to 2D.
    real(8), dimension(1:nsim, 0:lmax*(lmax+2)), device :: alm
    real(8), dimension(0:nring-1,-2:lmax      ), device :: plm
    
    ! scalar variables
    integer(4), intent(in)  :: m, type
    real(8),    intent(in)  :: beta
    integer                 :: stat
    real(8)  :: a1p, a2p
    ! compute FFT-matrix from alm, using dgemm provided by MKL
    ! IMPORTANT: get inputs from alm directly, no need for manual copying.

    stat = cublasSetStream(handle, stream)
    
    if (h_swap(1, type).eq.0) then
        stat =  cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, nsim, rlen, n, h_a1(1,type), alm(:,p1), nsim, plm(i0:i1,m), nring, beta, mat_cache1(:,i0), nsim)
        if (m.gt.0) &
             stat = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, nsim, rlen, n, h_a2(1,type), alm(:,p2), nsim, plm(i0:i1,m), nring, beta, mat_cache2(:,i0), nsim)
    else
        if (m.eq.0) then
            stat =  cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, nsim, rlen, n, h_a1(1,type), alm(:,p1), nsim, plm(i0:i1,m), nring, beta, mat_cache1(:,i0), nsim)           
        else
            stat =  cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, nsim, rlen, n, h_a1(1,type), alm(:,p1), nsim, plm(i0:i1,m), nring, beta, mat_cache2(:,i0), nsim)
            stat =  cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, nsim, rlen, n, h_a2(1,type), alm(:,p2), nsim, plm(i0:i1,m), nring, beta, mat_cache1(:,i0), nsim)
        endif
    endif
     
end subroutine

#else
SUBROUTINE alm2mat( alm, plm, m, type, beta)
    use sht_data_module, only: nside, nsim, npix, nring, lmax, mat_cache1, mat_cache2, &
                               i0, i1, p1, p2, n, rlen, swap, a1, a2
    implicit none

    ! Main input/output, note that alm is flattened from 3D to 2D.
    real(8), dimension(1:nsim, 0:lmax*(lmax+2)) :: alm
    real(8), dimension(0:nring-1,-2:lmax      ) :: plm

    ! scalar variables
    integer(4), intent(in)  :: m, type
    real(8),    intent(in)  :: beta

    ! compute FFT-matrix from alm, using dgemm provided by MKL
    ! IMPORTANT: get inputs from alm directly, no need for manual copying.
    if (swap(1, type).eq.0) then
        call dgemm("N", "T", nsim, rlen, n, a1(1,type), alm(:,p1), nsim, plm(i0:i1,m), nring, beta, mat_cache1(:,i0), nsim)
        if (m.gt.0) &
        call dgemm("N", "T", nsim, rlen, n, a2(1,type), alm(:,p2), nsim, plm(i0:i1,m), nring, beta, mat_cache2(:,i0), nsim)
    else
        if (m.eq.0) then
            call dgemm("N", "T", nsim, rlen, n, a1(1,type), alm(:,p1), nsim, plm(i0:i1,m), nring, beta, mat_cache1(:,i0), nsim)                
        else
            call dgemm("N", "T", nsim, rlen, n, a1(1,type), alm(:,p1), nsim, plm(i0:i1,m), nring, beta, mat_cache2(:,i0), nsim)
            call dgemm("N", "T", nsim, rlen, n, a2(1,type), alm(:,p2), nsim, plm(i0:i1,m), nring, beta, mat_cache1(:,i0), nsim)
        endif
    endif
end subroutine
#endif










!--------------------------------------------
! accessories
!--------------------------------------------
subroutine set_locations(m)
    use sht_data_module, only: i0, i1, n, rlen, p1, p2, ring_tab, lmax, nring
    implicit none
    integer(4), intent(in)  :: m

    n    = lmax+1-m     ! number of ells
    i0   = ring_tab(m)  ! For the ring optimization: start position
    rlen = nring-i0*2   ! For the ring optimization: length of the shortened ring
    i1   = i0+rlen-1    ! For the ring optimization: end position
    p1   = (lmax+2)*m   ! start position of the real-part in the flattened alm
    p2   = (lmax+1)*n   ! start position of the imaginary-part in the flattened alm

#ifdef GPU
    !$acc kernels
    n    = lmax+1-m     ! number of ells
    i0   = ring_tab(m)  ! For the ring optimization: start position
    rlen = nring-i0*2   ! For the ring optimization: length of the shortened ring
    i1   = i0+rlen-1    ! For the ring optimization: end position
    p1   = (lmax+2)*m   ! start position of the real-part in the flattened alm
    p2   = (lmax+1)*n   ! start position of the imaginary-part in the flattened alm
    !$acc end kernels
#endif
end subroutine set_locations



subroutine apply_mask_mini(map, mask, npix, nsim)
    implicit none
    real(8), dimension(0:npix-1, 1:nsim)    :: map
    real(8), dimension(0:npix-1)            :: mask
    integer(4)  :: npix, nsim, i

    do i=1,nsim
        map(:,i) = map(:,i)*mask
    enddo

end



subroutine copy_with_mask_mini(map1, map2, mask, npix, nsim)
    implicit none
    real(8), dimension(0:npix-1, 1:nsim)    :: map1, map2
    real(8), dimension(0:npix-1)            :: mask
    integer(4)  :: npix, nsim, i

    do i=1,nsim
        map2(:,i) = map1(:,i)*mask
    enddo

end



! simple program to compute the residual of linear fit
subroutine lin_resi_mini(x, y, resi, n, nsim, list, nlist, a0_arr, a1_arr)
    implicit none
    integer(4)  :: n, nsim, i, j, nlist, list(1:nlist)
    real(8)     :: x(0:n-1,1:nsim), y(0:n-1,1:nsim), resi(0:n-1,1:nsim), a0_arr(1:nsim), a1_arr(1:nsim)
    real(8)     :: mx, my, cxx, cxy

    do i=1,nsim
        mx  = sum(x(list,i))/nlist
        my  = sum(y(list,i))/nlist

        cxx = sum( (x(list,i)-mx)**2 )
        cxy = sum( (y(list,i)-my)*(x(list,i)-mx) )
        
        a1_arr(i)  = cxy/cxx 
        a0_arr(i)  = my - mx*a1_arr(i)

        resi(list,i)  = y(list,i) - a0_arr(i) - a1_arr(i)*x(list,i)
    enddo

    return 
END



! simple program to compute the residual of linear fit without mask
subroutine lin_resi_mini_nomask(x, y, n, nsim, a1_arr)
    implicit none
    integer(4)  :: n, nsim, i, j
    real(8)     :: x(0:n-1,1:nsim), y(0:n-1,1:nsim)
    real(8)     :: mx, my, cxx, cxy, a0, a1, a1_arr(1:nsim)

    do i=1,nsim
        mx  = sum(x(:,i))/n
        my  = sum(y(:,i))/n

        cxx = sum( (x(:,i)-mx)*(x(:,i)-mx) )
        cxy = sum( (y(:,i)-my)*(x(:,i)-mx) )
        a1  = cxy/cxx 
        a1_arr(i) = a1
    enddo

    return 
END
