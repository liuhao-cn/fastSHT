!
! contains: 
! cl2alm_TEB
! cl2alm_T
!

include 'mkl_vsl.f90'

SUBROUTINE cl2alm_TEB(argc, argv)         ! interface to IDL  
    !DEC$ ATTRIBUTES DLLEXPORT::cl2alm_TEB
    ! 
    INTEGER*8   argc, argv(*)                 ! Argc and Argv are L64 integers    

    CALL cl2alm_TEB1( %VAL(argv(1)),  %VAL(argv(2)),  %VAL(argv(3)), &
                      %VAL(argv(4)),  %VAL(argv(5)),  %VAL(argv(6)), &
                      %VAL(argv(7)),  %VAL(argv(8)) )

    RETURN     
END

! generate simulated alm_TEB from input cl, requires all 6 components: 
! TT, EE, BB, TE, TB, EB
!  1   2   3   4   5   6
! Empty components can be set to zero in the main program
SUBROUTINE cl2alm_TEB1( cl, lmin, lmax, nsim, almT, almE, almB, seed )
    USE MKL_VSL_TYPE
    USE MKL_VSL
    implicit none

    real(8)  :: cl(0:lmax, 1:6)
    real(8)  :: almT(1:nsim, 0:lmax, 0:lmax)
    real(8)  :: almE(1:nsim, 0:lmax, 0:lmax)
    real(8)  :: almB(1:nsim, 0:lmax, 0:lmax)
    ! other variables
    real(8), allocatable    :: a(:), b(:), c(:), u(:), v(:), w(:)
    real(8)                 :: fac
    integer(4)              :: nsim, lmin, lmax, l, m, i
    integer(4)              :: errcode, brng, method, seed
    TYPE (VSL_STREAM_STATE) :: stream

    allocate( a(0:lmax) ); a=0
    allocate( b(0:lmax) ); b=0
    allocate( c(0:lmax) ); c=0
    allocate( u(0:lmax) ); u=0 
    allocate( v(0:lmax) ); v=0
    allocate( w(0:lmax) ); w=0

    ! compute the matrix elements for each ell. They are used for linear
    ! combination to give the correct cross spectra
    do l=lmin, lmax
        a(l) = sqrt( cl(l,1) )
        b(l) = cl(l,4) / a(l)
        c(l) = cl(l,5) / a(l)
        u(l) = sqrt( cl(l,2) - b(l)**2 )
        v(l) = (cl(l,1)*cl(l,6) - cl(l,4)*cl(l,5))/(u(l)*cl(l,1))
        w(l) = sqrt( cl(l,3) - cl(l,5)/cl(l,1)/u(l)**2*(cl(l,2)*cl(l,5)-cl(l,4)*cl(l,6)) - v(l)*cl(l,6)/u(l) )
    enddo

    ! initialize the MKL random number generator
    brng = VSL_BRNG_MCG31
    method = VSL_RNG_METHOD_GAUSSIAN_ICDF
    errcode = vslnewstream(stream, brng,  seed)

    ! fill up the alms with initial Gaussian random numbers
    do m=0,lmax
        if (m.eq.0) then
            fac = 1
        else
            fac = sqrt(2.d0)/2
        endif
        errcode = vdrnggaussian( method, stream, nsim*(lmax+1), almT(:,:,m), 0.d0, fac)
        errcode = vdrnggaussian( method, stream, nsim*(lmax+1), almE(:,:,m), 0.d0, fac)
        errcode = vdrnggaussian( method, stream, nsim*(lmax+1), almB(:,:,m), 0.d0, fac)    
    enddo

    ! make correct alms by linear combinations. IMPORTATN: do it for B first, and
    ! then E, lastly T, so the values will not be updated too early
    do m=0,lmax
        do l=max(m, lmin),lmax
            almB(:,l,m) = almT(:,l,m)*c(l) + almE(:,l,m)*v(l) + almB(:,l,m)*w(l)
            almE(:,l,m) = almT(:,l,m)*b(l) + almE(:,l,m)*u(l)
            almT(:,l,m) = almT(:,l,m)*a(l)
            if (m.gt.0) then
                almB(:,l-m,lmax+1-m) = almT(:,l-m,lmax+1-m)*c(l) + almE(:,l-m,lmax+1-m)*v(l) + almB(:,l-m,lmax+1-m)*w(l)
                almE(:,l-m,lmax+1-m) = almT(:,l-m,lmax+1-m)*b(l) + almE(:,l-m,lmax+1-m)*u(l)
                almT(:,l-m,lmax+1-m) = almT(:,l-m,lmax+1-m)*a(l)
            endif
        enddo
    enddo

    ! zero the elements below lmin
    if (lmin.gt.0) then
        do m=0,lmin-1
            do l=m,lmin-1
                almB(:,l,m) = 0
                almE(:,l,m) = 0
                almT(:,l,m) = 0
                if (m.gt.0) then
                    almB(:,l-m,lmax+1-m) = 0
                    almE(:,l-m,lmax+1-m) = 0
                    almT(:,l-m,lmax+1-m) = 0
                endif
            enddo
        enddo
    endif

    return 
end






SUBROUTINE cl2alm_T(argc, argv)           ! interface to IDL  
    !DEC$ ATTRIBUTES DLLEXPORT::cl2alm_T
    ! 
    INTEGER*8   argc, argv(*)                 ! Argc and Argv are L64 integers    

    CALL cl2alm_T1( %VAL(argv(1)),  %VAL(argv(2)),  %VAL(argv(3)), &
                    %VAL(argv(4)),  %VAL(argv(5)),  %VAL(argv(6)) )

    RETURN     
END    

! generate simulated alm_TEB from input cl, only T, so do cl. 
SUBROUTINE cl2alm_T1( cl, lmin, lmax, nsim, almT, seed )
    USE MKL_VSL_TYPE
    USE MKL_VSL
    implicit none

    real(8)  :: cl(0:lmax)
    real(8)  :: almT(1:nsim, 0:lmax, 0:lmax)
    ! other variables
    real(8)                 :: fac
    integer(4)              :: nsim, lmin, lmax, l, m, i
    integer(4)              :: errcode, brng, method, seed
    TYPE (VSL_STREAM_STATE) :: stream

    ! initialize the MKL random number generator
    brng = VSL_BRNG_MCG31
    method = VSL_RNG_METHOD_GAUSSIAN_ICDF
    errcode = vslnewstream(stream, brng,  seed)

    ! for m=0
    do l=lmin,lmax
        errcode = vdrnggaussian( method, stream, nsim, almT(:,l,0), 0.d0, sqrt(cl(l)) )
    enddo

    ! for m>0
    do m=1,lmax
        do l=max(m, lmin),lmax
            errcode = vdrnggaussian( method, stream, nsim, almT(:,  l,       m), 0.d0, sqrt(cl(l)/2) )  
            errcode = vdrnggaussian( method, stream, nsim, almT(:,l-m,lmax+1-m), 0.d0, sqrt(cl(l)/2) )
        enddo
    enddo

    ! zero the elements below lmin
    if (lmin.gt.0) then
        do m=0,lmin-1
            do l=m,lmin-1
                almT(:,l,m) = 0
                if (m.gt.0) then
                    almT(:,l-m,lmax+1-m) = 0
                endif
            enddo
        enddo
    endif

    return 
end