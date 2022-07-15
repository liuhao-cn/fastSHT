!
! contains: 
! convert_alm
! convert_alm_healpy
! convert_alm_inv
! convert_alm_healpy_inv
!

! convert the new alm to the HEALPix fortran format
SUBROUTINE convert_alm1( alm_new, alm_trad, nsim, lmax)
    implicit none
    ! Main input/output, alm is saved as real-value square matrix with the new scheme (see fig.1)
    real(8)    :: alm_new(1:nsim, 0:lmax, 0:lmax)
    complex(8) :: alm_trad(0:lmax, 0:lmax, 1:nsim)
    ! other variables
    integer(4) :: nsim, lmax, m, i 

    do i=1,nsim
        alm_trad(:,0,i) = dcmplx(alm_new(i,:,0), alm_new(i,:,0)*0)
    enddo
     
    do m=1,lmax
    	do i=1,nsim
            alm_trad(m:lmax,m,i) = dcmplx(alm_new(i,m:lmax,m), alm_new(i,0:lmax-m,lmax+1-m))
        enddo
    enddo

    return 
end



! convert the new format alm to match the healpy format
SUBROUTINE convert_alm_healpy1( alm_new, alm_hpy, nsim, lmax)
    implicit none

    real(8)    :: alm_new(1:nsim, 0:lmax, 0:lmax)
    real(8)    :: alm_hpy(1:2, 1:(lmax+1)*(lmax+2)/2, 1:nsim)
    integer(4) :: nsim, lmax, l, m, i
    integer(8) :: pos 

    pos = 1
    do m=0,lmax
    	do l=m,lmax
            alm_hpy(1,pos,:) = alm_new(:,l,m)
            if (m.gt.0) then
                alm_hpy(2,pos,:) = alm_new(:,l-m,lmax+1-m)
            else
                alm_hpy(2,pos,:) = 0
            endif
            pos = pos + 1
    	enddo
    enddo

    return 
end



! inverse conversion of convert_alm
SUBROUTINE convert_alm_inv1( alm_trad, alm_new, nsim, lmax)
    implicit none
    ! Main input/output, alm is saved as real-value square matrix with the new scheme (see fig.1)
    real(8)    :: alm_new(1:nsim, 0:lmax, 0:lmax)
    complex(8) :: alm_trad(0:lmax, 0:lmax, 1:nsim)
    ! other variables
    integer(4) :: nsim, lmax, m, i 

    do i=1,nsim
        alm_new(i,:,0) = real(alm_trad(:,0,i))
    enddo
     
    do m=1,lmax
        do i=1,nsim
            alm_new(i,m:lmax,m)          = real ( alm_trad(m:lmax,m,i) )
            alm_new(i,0:lmax-m,lmax+1-m) = aimag( alm_trad(m:lmax,m,i) )
        enddo
    enddo

    return 
end




! inverse conversion of convert_alm_healpy
SUBROUTINE convert_alm_healpy_inv1( alm_hpy, alm_new, nsim, lmax)
    implicit none

    real(8)    :: alm_new(1:nsim, 0:lmax, 0:lmax)
    real(8)    :: alm_hpy(1:2, 1:(lmax+1)*(lmax+2)/2, 1:nsim)
    integer(4) :: nsim, lmax, l, m, i
    integer(8) :: pos 

    pos = 1
    do m=0,lmax
        do l=m,lmax
            alm_new(:,l,m) = alm_hpy(1,pos,:)
            if (m.gt.0) then
                alm_new(:,l-m,lmax+1-m) = alm_hpy(2,pos,:)
            endif
            pos = pos + 1
        enddo
    enddo

    return 
end
