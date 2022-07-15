SUBROUTINE norm_alm(argc, argv)           ! interface to IDL  
    !DEC$ ATTRIBUTES DLLEXPORT::norm_alm
    ! 
    INTEGER*8   argc, argv(*)             ! Argc and Argv are L64 integers    

    CALL norm_alm1( %VAL(argv(1)),  %VAL(argv(2)),  %VAL(argv(3)), &
                    %VAL(argv(4)),  %VAL(argv(5)),  %VAL(argv(6)), &
                    %VAL(argv(7)) )

    RETURN     
END    

! normalize the alm by a cl-like factor & shift: alm1 = alm*fac + shift
SUBROUTINE norm_alm1( alm_in, alm_out, fac, shift, lmax, nsim, m0 )
    implicit none
    real(8)     :: alm_in(1:nsim, 0:lmax, 0:lmax), alm_out(1:nsim, 0:lmax, 0:lmax)
    real(8)     :: fac(0:lmax), shift(0:lmax), m0
    integer(4)  :: nsim, lmax, l, m, i

    ! for m=0
    do l=0,lmax
        alm_out(:,l,0) = alm_in(:,l,0)*fac(l)*m0 + shift(l)
    enddo

    ! for m>0
    do m=1,lmax
        do l=m,lmax
            alm_out(:,  l,       m) = alm_in(:,  l,       m)*fac(l) + shift(l)
            alm_out(:,l-m,lmax+1-m) = alm_in(:,l-m,lmax+1-m)*fac(l) + shift(l)
        enddo
    enddo

    return 
end