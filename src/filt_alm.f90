!
! contains: 
! filt_alm
!
SUBROUTINE filt_alm(argc, argv)               ! interface to IDL  
    !DEC$ ATTRIBUTES DLLEXPORT::filt_alm
    INTEGER*8   argc, argv(*)                 ! Argc and Argv are L64 integers    

    call filt_alm1( %VAL(argv( 1)), %VAL(argv( 2)), %VAL(argv( 3)), &
                    %VAL(argv( 4)) )

    RETURN     
END

! filter the alm with a window function
SUBROUTINE filt_alm1( alm, win, nsim, lmax)
    implicit none
    ! Main input/output, alm is saved as real-value square matrix with the new scheme (see fig.1)
    real(8)    :: alm(1:nsim, 0:lmax, 0:lmax), win(0:lmax)
    ! other variables
    integer(4) :: nsim, lmax, ell, m

    ! for m=0
    do ell=0,lmax
        alm(:,ell,0) = alm(:,ell,0) * win(ell)
    enddo
    
    ! for m>0
    do m=1,lmax
        do ell=m,lmax
            alm(:,ell,m) = alm(:,ell,m) * win(ell)
            alm(:,ell-m,lmax+1-m) = alm(:,ell-m,lmax+1-m) * win(ell)
        enddo
    enddo

    return 
end