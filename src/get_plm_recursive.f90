!
! contains:
!
! get_plm_recursive: recursively compute plm

! pol=0: only temperature (plm0) 
!
! pol=1: only polarization (plm1 and plm2)
!
! note that m is fixed inside this program (as intent(in))
! 
! Part of the code is originated from the HEALPix program plm_gen.
!



subroutine get_plm_recursive(pol, m)
    use sht_data_module
    implicit none

    ! number-like constants for iterations
    integer(4), intent(in)    :: pol, m
    integer(4)	              :: l, i, m1
    real(8)                   :: fl

    !---------------------------------
    ! start recursion
    !---------------------------------

    !$acc kernels
    a(m) = 0
    b(m) = -1
    !$acc loop
    do l=m+1,lmax
       a(l)  = sqrt(  c3(l)/(c4(l)-m**2) )
       b(l)  = sqrt( (c1(l)-m**2)/c2(l)  ) * a(l)
    enddo


    ! use different scaling schemes for different resolutions
    if (nside.le.512) then
       ! for nside<=512, there is no scaling, only ring optimization
       !$acc loop
        do i=i0,nhr
            plm0(i,m-1) = 0
            plm0(i,m-2) = sqrt( diag_fac(m)/4/pi )*(-ss(i))**m
        enddo

        !$acc loop gang
        do i=i0,nhr
           !$acc loop vector seq
           do l=m,lmax
              plm0(i,l) = a(l)*cc(i)*plm0(i,l-1) - b(l)*plm0(i,l-2)
           enddo
        enddo
        ! plm0(i0:nhr,m) = sqrt( diag_fac(m)/4/pi )*(-ss(i0:nhr))**m

        ! ! continue to the non-diagonal plms (intensity)
        ! if (m.lt.lmax) then        
        !     ! the l=m+1 element
        !     plm0(i0:nhr,m+1) = a(m+1)*cc(i0:nhr)*plm0(i0:nhr,m)

        !     ! other elements
        !     do l=m+2, lmax
        !         plm0(i0:nhr,l) = a(l)*cc(i0:nhr)*plm0(i0:nhr,l-1) - b(l)*plm0(i0:nhr,l-2)
        !     enddo
        ! endif
    else
        ! for nside>512, start from pre-computed safe positions to avoid scaling.
        ! note that the settings ensure that "safe_pos" is less than lmax, so
        ! there is no need to use "if" in the core loops.
        if (m.eq.lmax) then
            plm0(i0:nhr,lmax) = sqrt( diag_fac(m)/(4*pi) )*(-ss(i0:nhr))**lmax
         else
            !$acc loop
            do i=i0,nhr 
                ! if(pos(i,m)-1 .ge. m) plm0(i,m:pos(i,m)-1) = 0
                plm0(i,  pos(i,m)  ) = val1(i,m)
                plm0(i,  pos(i,m)+1) = val2(i,m)
            enddo

            !$acc loop gang
            do i=i0,nhr
               !$acc loop vector seq
               do l=m, lmax
                  if(l.ge.pos(i,m)+2) plm0(i,l) = a(l)*cc(i)*plm0(i,l-1) - b(l)*plm0(i,l-2)
                enddo
            enddo
        endif
    endif


    !$acc loop independent
    do i=nside*2,nring-1-i0                
       do l=m,lmax
            plm0(i,l) = plm0(nring-1-i,l)*(-1)**(l+m)
        enddo
    enddo



    ! for the pol-case
    if (pol.eq.1) then
        ! compute lam_fact and c_w
        m1  = max(1,m)  
        lam_fact(0:m1) = 0
        if (m1.lt.lmax) forall(l=m1+1:lmax) lam_fact(l) = 2*sqrt( dble(2*l+1)/(2*l-1)*(l**2-m**2) )
        c_w(i0:nhr) = m/ss(i0:nhr)**2

        ! compute the diagonal elements of pol-plm
        plm1(i0:nhr,m) = normal_m(m)*normal_l(m)*plm0(i0:nhr,m)*(0.5d0-1.d0/ss(i0:nhr)**2)
        plm2(i0:nhr,m) = normal_m(m)*normal_l(m)*plm0(i0:nhr,m)*cc(i0:nhr)/ss(i0:nhr)**2

        ! compute the non-diagonal terms of pol-plm
        if (m.lt.lmax) then
           !$acc loop
           do i=i0,nhr
              do l=m+1, lmax
                 !fl = l
                    plm1(i,l) = ( plm0(i,l-1)*lam_fact(l)*b_w(i) - plm0(i,l)*(two_on_s2(i)*(l-m**2) +l*(l-1.d0)) )*normal_l(l)
                    plm2(i,l) = ( plm0(i,l-1)*lam_fact(l)        - plm0(i,l)*(2.d0*cc(i)*(l-1.d0)) ) *c_w(i)*normal_l(l)
                enddo
            enddo
        endif

        ! extend to south hemisphere
        !$acc loop independent
        do l=m,lmax
            do i=nside*2,nring-1-i0
                plm1(i,l) = plm1(nring-1-i,l)*(-1)**(l+m)
                plm2(i,l) = plm2(nring-1-i,l)*(-1)**(l+m+1)
            enddo
        enddo
    endif
    !$acc end kernels
end subroutine











!---------------------------
! obsolete code
!---------------------------!
! ! compute the diagonal plms (l=m)
!     if (nside.lt.1024) then
!         ! for nside=1~512: no scaling
!         plm0(i0:nhr,m) = sqrt( diag_fac/4/pi )*(-ss(i0:nhr))**m
!     else
!         ! for nside=1024: use simple scaling : apply 50% of sin(theta)^m at the
!         ! beginning and end respectively, and do not scale the "safe" rings. 
!         ! borrow the memory of c_w
!         rth = 1023
!         if (i0.le.rth) then
!             c_w(i0:rth) = ss(i0:rth)**(m-m/2)
!             c_w(rth+1:nhr) = 1
!             ! make the scaled diagonal term
!             plm0(i0:rth, m)   = sqrt( diag_fac/4/pi ) * (-1)**m * ss(i0:rth )**(m/2)
!             plm0(rth+1:nhr,m) = sqrt( diag_fac/4/pi ) * (-1)**m * ss(rth+1:nhr)**m
!         else
!             plm0(i0:nhr,m) = sqrt( diag_fac/4/pi )*(-ss(i0:nhr))**m
!         endif
!     endif
!
!     ! for nside=1~1024: continue to the non-diagonal plms (intensity)
!     if (m.lt.lmax .and. nside.le.1024) then        
!         ! the l=m+1 element
!         plm0(i0:nhr,m+1) = a(m+1)*cc(i0:nhr)*plm0(i0:nhr,m)
!
!         ! other elements
!         do l=m+2, lmax
!             plm0(i0:nhr,l) = a(l)*cc(i0:nhr)*plm0(i0:nhr,l-1) - b(l)*plm0(i0:nhr,l-2)
!         enddo
!     endif
!
!     ! for nside=1024: apply another 50% of sin(theta)^m 
!     if (nside.eq.1024 .and. i0.le.rth) forall(l=m:lmax) plm0(i0:rth,l) = plm0(i0:rth,l)*c_w(i0:rth)
