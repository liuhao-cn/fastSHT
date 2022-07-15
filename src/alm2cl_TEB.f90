!
! contains: 
! alm2cl_TEB
! alm2cl_T
!

SUBROUTINE alm2cl_TEB(argc, argv)         ! interface to IDL  
    !DEC$ ATTRIBUTES DLLEXPORT::alm2cl_TEB
    ! 
    INTEGER*8   argc, argv(*)                 ! Argc and Argv are L64 integers    

    CALL alm2cl_TEB1( %VAL(argv(1)),  %VAL(argv(2)),  %VAL(argv(3)), &
                      %VAL(argv(4)),  %VAL(argv(5)),  %VAL(argv(6)), &
                      %VAL(argv(7)) )

    RETURN     
END    

! compute cl from the new format alm, all 6 components
! TT, EE, BB, TE, TB, EB
!  1   2   3   4   5   6
SUBROUTINE alm2cl_TEB1( cl, lmin, lmax, nsim, almT, almE, almB )
    implicit none

    real(8)  :: cl(1:nsim, 0:lmax, 1:6)
    real(8)  :: almT(1:nsim, 0:lmax, 0:lmax)
    real(8)  :: almE(1:nsim, 0:lmax, 0:lmax)
    real(8)  :: almB(1:nsim, 0:lmax, 0:lmax)
    ! other variables
    real(8)                 :: fac
    integer(4)              :: nsim, lmin, lmax, l, m, l1

    ! compute cl from alm
    do m=0,lmax
        if (m.eq.0) then
            fac = 1
        else
            fac = 2
        endif
        l1 = max(m, lmin)

        do l=l1,lmax
            ! real-part
            cl(:,l,1) = cl(:,l,1) + almT(:,l,m)**2*fac/(2*l+1)
            cl(:,l,2) = cl(:,l,2) + almE(:,l,m)**2*fac/(2*l+1)
            cl(:,l,3) = cl(:,l,3) + almB(:,l,m)**2*fac/(2*l+1)
            cl(:,l,4) = cl(:,l,4) + almT(:,l,m)*almE(:,l,m)*fac/(2*l+1)
            cl(:,l,5) = cl(:,l,5) + almT(:,l,m)*almB(:,l,m)*fac/(2*l+1)
            cl(:,l,6) = cl(:,l,6) + almE(:,l,m)*almB(:,l,m)*fac/(2*l+1)
            ! imaginary-part
            if (m.gt.0) then
                cl(:,l,1) = cl(:,l,1) + almT(:,l-m,lmax+1-m)**2*fac/(2*l+1)
                cl(:,l,2) = cl(:,l,2) + almE(:,l-m,lmax+1-m)**2*fac/(2*l+1)
                cl(:,l,3) = cl(:,l,3) + almB(:,l-m,lmax+1-m)**2*fac/(2*l+1)
                cl(:,l,4) = cl(:,l,4) + almT(:,l-m,lmax+1-m)*almE(:,l-m,lmax+1-m)*fac/(2*l+1)
                cl(:,l,5) = cl(:,l,5) + almT(:,l-m,lmax+1-m)*almB(:,l-m,lmax+1-m)*fac/(2*l+1)
                cl(:,l,6) = cl(:,l,6) + almE(:,l-m,lmax+1-m)*almB(:,l-m,lmax+1-m)*fac/(2*l+1)
            endif
        enddo
    enddo

    return 
end




SUBROUTINE alm2cl_T(argc, argv)           ! interface to IDL  
    !DEC$ ATTRIBUTES DLLEXPORT::alm2cl_T
    ! 
    INTEGER*8   argc, argv(*)                 ! Argc and Argv are L64 integers    

    CALL alm2cl_T1( %VAL(argv(1)),  %VAL(argv(2)),  %VAL(argv(3)), &
                    %VAL(argv(4)),  %VAL(argv(5)) )

    RETURN     
END    


! compute cl from alm, only for T
SUBROUTINE alm2cl_T1( cl, lmin, lmax, nsim, almT)
    implicit none

    real(8)  :: cl(1:nsim, 0:lmax)
    real(8)  :: almT(1:nsim, 0:lmax, 0:lmax)
    real(8)  :: almE(1:nsim, 0:lmax, 0:lmax)
    real(8)  :: almB(1:nsim, 0:lmax, 0:lmax)
    ! other variables
    real(8)                 :: fac
    integer(4)              :: nsim, lmin, lmax, l, m, l1

    ! compute cl from alm
    do m=0,lmax
        if (m.eq.0) then
            fac = 1
        else
            fac = 2
        endif
        l1 = max(m, lmin)

        do l=l1,lmax
            ! real-part
            cl(:,l) = cl(:,l) + almT(:,l,m)**2*fac/(2*l+1)
            ! imaginary-part
            if (m.gt.0) cl(:,l) = cl(:,l) + almT(:,l-m,lmax+1-m)**2*fac/(2*l+1)
        enddo
    enddo

    return 
end
