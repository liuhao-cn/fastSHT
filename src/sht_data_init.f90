! allocate the variables in "sht_data_module", actual code
subroutine sht_data_alloc(params)
  use sht_data_module
#ifdef GPU
  use cublas_v2
#endif
    implicit none
    integer(4) :: params(1:5), nbuff, stat, m, i, k

    nside    = params(1)
    lmax     = params(2)
    nring    = params(3)
    nsim     = params(4)
    nbuff    = params(5)
    is_pol      = params(6)
    
    npix     = 12*int(nside)**2
    pa       = 4*pi/npix
    nhr      = 2*nside-1
    i0       = 0
    i1       = 0
    p1       = 0
    p2       = 0
    n        = 0
    rlen     = 0

#ifdef GPU
    !$acc kernels
    nside    = params(1)
    lmax     = params(2)
    nring    = params(3)
    nsim     = params(4)
    nbuff    = params(5)

    npix     = 12*int(nside)**2
    pa       = 4*pi/npix
    nhr      = 2*nside-1
    i0       = 0
    i1       = 0
    p1       = 0
    p2       = 0
    n        = 0
    rlen     = 0
    !$acc end kernels
    stat = cublasCreate(handle)
    stat = cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_DEVICE)

#endif
    
    ! ! ring table and safe positions
    if ( allocated(ring_tab ) ) deallocate(ring_tab)
    if ( allocated(val1     ) ) deallocate(val1    )
    if ( allocated(val2     ) ) deallocate(val2    )
    if ( allocated(pos      ) ) deallocate(pos     )

    allocate( ring_tab(            0:lmax) ); ring_tab = 0;
    allocate( val1    (0:2*nside-1,0:lmax) ); val1 = 0;
    allocate( val2    (0:2*nside-1,0:lmax) ); val2 = 0;
    allocate( pos     (0:2*nside-1,0:lmax) ); pos = 0;
    
    ! basic ring information
    if (allocated(theta))  deallocate(theta)
    if (allocated(phi0 ))  deallocate(phi0 )

    allocate( theta(0:nside*2-1) ); theta = 0;
    allocate( phi0 (0:nside*2-1) ); phi0 = 0;

    ! initialization for the plm-constants
    if (allocated(cc       )) deallocate(cc )
    if (allocated(ss       )) deallocate(ss )
    if (allocated(b_w      )) deallocate(b_w)
    if (allocated(c_w      )) deallocate(c_w)
    if (allocated(two_on_s2)) deallocate(two_on_s2)

    allocate( cc       (0:nside*2-1) ); cc = 0;
    allocate( ss       (0:nside*2-1) ); ss = 0;
    allocate( b_w      (0:nside*2-1) ); b_w = 0;
    allocate( c_w      (0:nside*2-1) ); c_w = 0;
    allocate( two_on_s2(0:nside*2-1) ); two_on_s2 = 0;
    
    if (allocated(c1      )) deallocate(c1)
    if (allocated(c2      )) deallocate(c2)
    if (allocated(c3      )) deallocate(c3)
    if (allocated(c4      )) deallocate(c4)
    if (allocated(a       )) deallocate(a)
    if (allocated(b       )) deallocate(b)
    if (allocated(normal_l)) deallocate(normal_l)
    if (allocated(normal_m)) deallocate(normal_m)
    if (allocated(lam_fact)) deallocate(lam_fact)
    if (allocated(diag_fac )) deallocate(diag_fac)

    allocate( c1      (0:lmax) ); c1 = 0;
    allocate( c2      (0:lmax) ); c2 = 0;
    allocate( c3      (0:lmax) ); c3 = 0;
    allocate( c4      (0:lmax) ); c4 = 0;
    allocate( a       (0:lmax) ); a = 0;
    allocate( b       (0:lmax) ); b = 0;
    allocate( normal_l(0:lmax) ); normal_l = 0; 
    allocate( normal_m(0:lmax) ); normal_m = 0;
    allocate( lam_fact(0:lmax) ); lam_fact = 0;
    allocate( diag_fac (0:lmax) ); diag_fac = 0;

    if (allocated(plm0)) deallocate(plm0)
    if (allocated(plm1)) deallocate(plm1)
    if (allocated(plm2)) deallocate(plm2)

    allocate( plm0(0:nring-1,-2:lmax) ); plm0 = 0;
    allocate( plm1(0:nring-1,-2:lmax) ); plm1 = 0;
    allocate( plm2(0:nring-1,-2:lmax) ); plm2 = 0;
    
    if (allocated(i1_arr)) deallocate(i1_arr)
    if (allocated(i2_arr)) deallocate(i2_arr)
    if (allocated(nfft_arr)) deallocate(nfft_arr)

    allocate( i1_arr(0:nring-1)   ); i1_arr = 0;
    allocate( i2_arr(0:nring-1)   ); i2_arr = 0;
    allocate( nfft_arr(0:nring-1) ); nfft_arr = 0;
    

    
    ! compute the FFT-mapping table and factors
    ! Forward
    if (allocated(tab )) deallocate(tab)
    if (allocated(fac1)) deallocate(fac1)
    if (allocated(fac2)) deallocate(fac2)
    if (allocated(fac3)) deallocate(fac3)
    if (allocated(fac4)) deallocate(fac4)

    allocate( tab (0:  nring-1,0:lmax) ); tab = 0;
    allocate( fac1(0:2*nside-1,0:lmax) ); fac1 = 0;
    allocate( fac2(0:2*nside-1,0:lmax) ); fac2 = 0;
    allocate( fac3(0:2*nside-1,0:lmax) ); fac3 = 0;
    allocate( fac4(0:2*nside-1,0:lmax) ); fac4 = 0;

    ! Backward
    if (allocated(tabI )) deallocate(tabI)
    if (allocated(fac1I)) deallocate(fac1I)
    if (allocated(fac2I)) deallocate(fac2I)
    if (allocated(fac3I)) deallocate(fac3I)
    if (allocated(fac4I)) deallocate(fac4I)

    allocate( tabI (0:nring-1,0:lmax)   ); tabI = 0;
    allocate( fac1I(0:2*nside-1,0:lmax) ); fac1I = 0;
    allocate( fac2I(0:2*nside-1,0:lmax) ); fac2I = 0;
    allocate( fac3I(0:2*nside-1,0:lmax) ); fac3I = 0;
    allocate( fac4I(0:2*nside-1,0:lmax) ); fac4I = 0;
    
    ! tables for dgemm factors
    if (allocated(swap)) deallocate(swap)
    if (allocated(a1  )) deallocate(a1)
    if (allocated(a2  )) deallocate(a2)

    allocate( swap(0:1,0:3) ); swap = 0;
    allocate( a1  (0:1,0:3) ); a1 = 0;
    allocate( a2  (0:1,0:3) ); a2 = 0;

    
    ! slices of the FFT-matrix
    if (allocated(mat_cache1)) deallocate(mat_cache1)
    if (allocated(mat_cache2)) deallocate(mat_cache2)

    allocate( mat_cache1(1:nsim, 0:nring-1) ); mat_cache1 = 0;
    allocate( mat_cache2(1:nsim, 0:nring-1) ); mat_cache2 = 0;

    !main FFT-buffs for SHT
    if (nbuff.ge.1) then
        if (allocated(buff1)) deallocate(buff1)
        allocate( buff1(1:nsim, 0:npix-1) );
#ifdef GPU
        if (allocated(d_buff1)) deallocate(d_buff1)
        allocate( d_buff1(1:nsim, 0:npix-1) );
#else
        !$omp parallel do
        do i=0, npix-1
           do k=1,nsim
              buff1(k, i) = 0
           enddo
        enddo
#endif
    endif

    if (nbuff.ge.2) then
        if (allocated(buff2)) deallocate(buff2)
#ifdef GPU
        if (allocated(d_buff2)) deallocate(d_buff2)
        allocate( d_buff2(1:nsim, 0:npix-1) );
#else
        allocate( buff2(1:nsim, 0:npix-1) );
        !$omp parallel do
        do i=0, npix-1
           do k=1,nsim
              buff2(k, i) = 0
           enddo
        enddo
#endif
    endif

    if (nbuff.ge.3) then
        if (allocated(buff3)) deallocate(buff3)
#ifdef GPU
        if (allocated(d_buff3)) deallocate(d_buff3)
        allocate( d_buff3(1:nsim, 0:npix-1) );
#else
        allocate( buff3(1:nsim, 0:npix-1) );
        !$omp parallel do
        do i=0, npix-1
           do k=1,nsim
              buff3(k, i) = 0
           enddo
        enddo
#endif
    endif

    if (nbuff.ge.4) then
       if (allocated(buff4)) deallocate(buff4)
#ifdef GPU
       if (allocated(d_buff4)) deallocate(d_buff4)
       allocate( d_buff4(1:nsim, 0:npix-1) );
#else
       allocate( buff4(1:nsim, 0:npix-1) );
       !$omp parallel do
       do i=0, npix-1
          do k=1,nsim
             buff4(k, i) = 0
           enddo
        enddo
#endif
    endif

#ifdef GPU
    if (allocated(d_alm)) deallocate(d_alm)
    if (allocated(h_a1)) deallocate(h_a1)
    if (allocated(h_a2)) deallocate(h_a2)
    if (allocated(h_swap)) deallocate(h_swap)
    allocate( d_alm(1:nsim, 0:lmax, 0:lmax) ); d_alm = 0;
    allocate( h_a1  (0:1,0:3) ); h_a1 = 0;
    allocate( h_a2  (0:1,0:3) ); h_a2 = 0;
    allocate( h_swap(0:1,0:3) ); h_swap = 0;
    if (allocated(h_i1_arr)) deallocate(h_i1_arr)
    if (allocated(h_i2_arr)) deallocate(h_i2_arr)
    if (allocated(h_nfft_arr)) deallocate(h_nfft_arr)
    allocate( h_i1_arr(0:nring-1)  ); h_i1_arr = 0;
    allocate( h_i2_arr(0:nring-1)  ); h_i2_arr = 0;
    allocate( h_nfft_arr(0:nring-1) ); h_nfft_arr = 0;
    
    if(is_pol > 0) then
       if (allocated(d_alm2)) deallocate(d_alm2)
       allocate( d_alm2(1:nsim, 0:lmax, 0:lmax) );    d_alm2 = 0;
    endif

#endif
    ! flags for set/accumulate the map-FFT
    if (allocated(acc_flag1)) deallocate(acc_flag1)
    allocate( acc_flag1(0:npix-1) ); acc_flag1 = 0;

    if (allocated(acc_flag2)) deallocate(acc_flag2)
    allocate( acc_flag2(0:npix-1) ); acc_flag2 = 0;

#ifdef GPU
    if (allocated(cu_streams)) then
       deallocate(cu_streams)
    endif
    allocate(cu_streams(0:max(100, 4*nside )))
    do m=0, 4*nside
        stat = cudaStreamCreate(cu_streams(m))
    enddo
#endif
    
end subroutine





subroutine sht_set_data(ring_tab_in, pos_in, val1_in, val2_in, theta_in, phi0_in)
  use sht_data_module
#ifdef GPU
  use cudafor
#endif
    implicit none

    integer(4), intent(in), dimension(0:lmax)             :: ring_tab_in   
    integer(4), intent(in), dimension(0:2*nside-1,0:lmax) :: pos_in
    real(8),    intent(in), dimension(0:2*nside-1,0:lmax) :: val1_in, val2_in
    real(8),    intent(in), dimension(0:2*nside-1)        :: theta_in, phi0_in

    integer stat

    ring_tab = ring_tab_in
    pos = pos_in
    val1 = val1_in
    val2 = val2_in
    theta = theta_in
    phi0 = phi0_in

    ! compute the ring information, this must be done before computing the
    ! mapping tables
    call get_ring_info()

    ! compute the FFT mapping table and factors
    call fft2mat_compute_table()

    call mat2fft_compute_table()

    ! compute the recursive constants that do not change with m
    call set_out_of_loop_constants()

    ! compute the factors for dgemm
    call set_gemm_params_for_fft2alm()

end subroutine







subroutine sht_data_free()
    use sht_data_module
    implicit none

    if ( allocated(cc       ))  deallocate(cc       )
    if ( allocated(ss       ))  deallocate(ss       )
    if ( allocated(c1       ))  deallocate(c1       )
    if ( allocated(c2       ))  deallocate(c2       )
    if ( allocated(c3       ))  deallocate(c3       )
    if ( allocated(c4       ))  deallocate(c4       )
    if ( allocated(a        ))  deallocate(a        )
    if ( allocated(b        ))  deallocate(b        )
    if ( allocated(b_w      ))  deallocate(b_w      )
    if ( allocated(c_w      ))  deallocate(c_w      )
    if ( allocated(two_on_s2))  deallocate(two_on_s2)
    if ( allocated(normal_l ))  deallocate(normal_l )
    if ( allocated(normal_m ))  deallocate(normal_m )
    if ( allocated(lam_fact ))  deallocate(lam_fact )
    if ( allocated(diag_fac ))  deallocate(diag_fac )
    if ( allocated(theta    ))  deallocate(theta    )
    if ( allocated(phi0     ))  deallocate(phi0     )
    if ( allocated(plm0     ))  deallocate(plm0     )
    if ( allocated(plm1     ))  deallocate(plm1     )
    if ( allocated(plm2     ))  deallocate(plm2     )
    if ( allocated(i1_arr   ))  deallocate(i1_arr   )
    if ( allocated(i2_arr   ))  deallocate(i2_arr   )
    if ( allocated(nfft_arr ))  deallocate(nfft_arr )
    if ( allocated(swap     ))  deallocate(swap     )
    if ( allocated(a1       ))  deallocate(a1       )
    if ( allocated(a2       ))  deallocate(a2       )
    if ( allocated(tab      ))  deallocate(tab      )  
    if ( allocated(fac1     ))  deallocate(fac1     )  
    if ( allocated(fac2     ))  deallocate(fac2     )  
    if ( allocated(fac3     ))  deallocate(fac3     )  
    if ( allocated(fac4     ))  deallocate(fac4     ) 
    if ( allocated(tabI     ))  deallocate(tabI     )
    if ( allocated(fac1I    ))  deallocate(fac1I    ) 
    if ( allocated(fac2I    ))  deallocate(fac2I    ) 
    if ( allocated(fac3I    ))  deallocate(fac3I    ) 
    if ( allocated(fac4I    ))  deallocate(fac4I    ) 
    if ( allocated(ring_tab ))  deallocate(ring_tab )
    if ( allocated(val1     ))  deallocate(val1     )
    if ( allocated(val2     ))  deallocate(val2     )
    if ( allocated(pos      ))  deallocate(pos      )
    if ( allocated(buff1    ))  deallocate(buff1    )    
    if ( allocated(buff2    ))  deallocate(buff2    )    
    if ( allocated(buff3    ))  deallocate(buff3    )    
    if ( allocated(buff4    ))  deallocate(buff4    )
    if ( allocated(acc_flag1))  deallocate(acc_flag1)
    if ( allocated(acc_flag2))  deallocate(acc_flag2)
    if ( allocated(mat_cache1)) deallocate(mat_cache1)
    if ( allocated(mat_cache2)) deallocate(mat_cache2)

end subroutine sht_data_free









! pre-compute the recursion constant that do not change with m
subroutine set_out_of_loop_constants( )
    use sht_data_module, only: nside, lmax, cc, ss, theta, two_on_s2, &
         normal_l, normal_m, c1, c2, c3, c4, b_w, diag_fac
        
    implicit none

    integer(4)             :: l, m

    !$acc kernels
    cc = cos(theta)
    ss = sin(theta)
    two_on_s2 = 2/ss**2


    ! normalization factors for polarization
    normal_l(0:1) = 0 
    forall(l=2:lmax) normal_l(l) = 1/sqrt( dble(l+2)*dble(l+1)*dble(l)*dble(l-1) )
    forall(m=0:lmax) normal_m(m) = 2*m*(1-m)
    forall(l=1:lmax) c1(l) = (l-1)**2
    forall(l=1:lmax) c2(l) = (2*l-3)*(2*l-1)
    forall(l=1:lmax) c3(l) = (2*l-1)*(2*l+1)
    forall(l=1:lmax) c4(l) = l**2
    c1(0)=0; c2(0)=0; c3(0)=0; c4(0)=0;

    ! b_w is one of the factors for polarization
    b_w = cc/ss**2

    ! the diagonal elements normalization factor
    diag_fac(0) = 1

    !$acc loop seq
    do m=1, lmax
        diag_fac(m) = diag_fac(m-1)*(2*m+1)/2/m
     enddo
     !$acc end kernels

end subroutine













! Parse the Perm format real-value map-FFT, load values into a rectangular 
! FFT-matrix for m in [m1, m2]. Note that map-FFT is already transposed
!
! pre-compute the table and factor here, so the main action can be faster.
!
subroutine fft2mat_compute_table()
    use sht_data_module, only: nside, npix, nring, lmax, phi0, i1_arr, nfft_arr, tab,  fac1,  fac2,  fac3,  fac4
    implicit none

    real(8)     :: cc, ss
    integer(4)  :: m, pos, i, i1, nfft

    ! posietion
    !$acc kernels
    !$acc loop private(i1, nfft, pos)
    do i=0,nring-1       
        i1 = i1_arr(i)
        nfft = nfft_arr(i)                        ! nfft is the actual FFT-size
        !!$acc loop private(pos)
        do m=0,lmax    
            pos  = mod(m, nfft)
            if (pos.eq.0) then                    ! zero FFT-frequency       
                tab(i,m) = i1
            else if (pos.lt.nfft/2) then          ! positive FFT-frequency
                tab(i,m) = i1+pos*2
            else if (pos.eq.nfft/2) then          ! Nyquist FFT-frequency
                tab(i,m) = i1+1
            else                                  ! negative FFT-frequency            
                tab(i,m) = i1+(nfft-pos)*2
            endif
        enddo
    enddo

    ! factor
    !$acc loop private(i1, nfft, pos, cc, ss)
    do i=0,nside*2-1
        i1 = i1_arr(i)
        nfft = nfft_arr(i)                        ! nfft is the actual FFT-size
        !!$acc loop private(pos, cc, ss)
        do m=0,lmax    
            pos  = mod(m, nfft)
            cc   = cos( -m*phi0(i) )
            ss   = sin( -m*phi0(i) )
            if (pos.eq.0) then                    ! zero FFT-frequency       
                fac1(i,m) = cc; fac2(i,m) = 0
                fac3(i,m) = ss; fac4(i,m) = 0
            else if (pos.lt.nfft/2) then          ! positive FFT-frequency
                fac1(i,m) = cc; fac2(i,m) = -ss
                fac3(i,m) = ss; fac4(i,m) = cc
            else if (pos.eq.nfft/2) then          ! Nyquist FFT-frequency
                fac1(i,m) = cc; fac2(i,m) = 0
                fac3(i,m) = ss; fac4(i,m) = 0
            else                                  ! negative FFT-frequency            
                fac1(i,m) = cc; fac2(i,m) = ss
                fac3(i,m) = ss; fac4(i,m) =-cc
            endif
        enddo
    enddo
    !$acc end kernels
end subroutine













! Compute the table and factors used to map the FFT-matrix to the map-FFT
! Detailed illustrations are in the old "mat2fft_run" subroutine.   
subroutine mat2fft_compute_table()
    use sht_data_module, only: npix, nring, lmax, phi0, i1_arr, nfft_arr, tabI, fac1I, fac2I, fac3I, fac4I
    implicit none

    real(8)     :: cc, ss
    integer(4)  :: m, pos, i, i1, nside, nfft, fac_m0

    !$acc kernels

    nside = nint( (npix/12)**0.5, 4 ) 

    ! position
    !$acc loop private(i1, nfft, pos)
    do i=0,nring-1       
        i1 = i1_arr(i)
        nfft = nfft_arr(i)                        
        do m=0,lmax
            pos  = mod( m, nfft)
            if (pos.eq.0) then                      
                ! zero frequency
                tabI(i,m) = i1     
            else if (pos.lt.nfft/2) then            
                ! positive frequency (already in the upper unit circle), +1 for imaginary part
                tabI(i,m) = i1+pos*2
            else if (pos.eq.nfft/2) then            
                ! Nyquist frequency, in the Perm format (even size), this is next to zero frequency
                tabI(i,m) = i1+1
            else                                    
                ! negative frequency, take conjugate and flip to the upper half of the unit circle
                tabI(i,m) = i1+(nfft-pos)*2
            endif             
        enddo
    enddo
    ! factor
    !$acc loop private(i1, nfft, pos, cc, ss, fac_m0)
    do i=0,2*nside-1      

        i1 = i1_arr(i)
        nfft = nfft_arr(i)                        

        do m=0,lmax
            if (m.eq.0) fac_m0 = 1
            if (m.ne.0) fac_m0 = 2
 
            pos  = mod( m, nfft)
            cc   = cos( m*phi0(i) )
            ss   = sin( m*phi0(i) )
            if (pos.eq.0) then                      
                ! zero frequency
                fac1I(i,m) = cc*fac_m0; fac2I(i,m) = -ss*fac_m0
                fac3I(i,m) = 0;               fac4I(i,m) = 0       
            else if (pos.lt.nfft/2) then            
                ! positive frequency (already in the upper unit circle), +1 for imaginary part
                fac1I(i,m) = cc;       fac2I(i,m) = -ss
                fac3I(i,m) = ss;       fac4I(i,m) =  cc
            else if (pos.eq.nfft/2) then            
                ! Nyquist frequency, in the Perm format (even size), this is next to zero frequency
                fac1I(i,m) = cc*fac_m0; fac2I(i,m) = -ss*fac_m0
                fac3I(i,m) = 0;              fac4I(i,m) =  0
            else                                    
                ! negative frequency, take conjugate and flip to the upper half of the unit circle
                fac1I(i,m) = cc;       fac2I(i,m) = -ss
                fac3I(i,m) =-ss;       fac4I(i,m) = -cc
            endif             
        enddo
    enddo
    !$acc end kernels
end subroutine











! Compute the constants that should be used in calling "dgemm"
! 
! "a1", "a2" are the constants to be computed, determined by "type"
!
! "swap" controls the swap of real/imaginary parts. 
!
! This helps "dgemm" to swallow as much job as possible.
! 
subroutine set_gemm_params_for_fft2alm( )
#ifdef GPU
  use sht_data_module, only: nside, pi, swap, a1, a2, pa, h_swap, h_a1, h_a2
#else
  use sht_data_module, only: nside, pi, swap, a1, a2, pa
#endif
    implicit none
    integer(4) :: inv, type

    !$acc kernels
    !$acc loop
    do inv=0,1
        do type=0,3
            if (inv.eq.0) then             ! settings for the forward transform
                select case (type)
                  case(0)                  ! accumulate alm directly, the simplest case, for intensity
                    a1(inv,type)= pa; a2(inv,type)= pa; swap(inv,type)=0 
                  case(1)                  ! accumulate -alm
                    a1(inv,type)=-pa; a2(inv,type)=-pa; swap(inv,type)=0 
                  case(2)                  ! accumulate i*alm
                    a1(inv,type)= pa; a2(inv,type)=-pa; swap(inv,type)=1 
                  case default             ! accumulate -i*alm
                    a1(inv,type)=-pa; a2(inv,type)= pa; swap(inv,type)=1 
                end select
            else                           ! settings for the backward transform
                select case (type)
                  case(0)                  ! use alm directly, the simplest case, for intensity
                    a1(inv, type)= 1; a2(inv, type)= 1; swap(inv, type)=0   
                  case(1)                  ! use -alm
                    a1(inv, type)=-1; a2(inv, type)=-1; swap(inv, type)=0   
                  case(2)                  ! use i*alm
                    a1(inv, type)= 1; a2(inv, type)=-1; swap(inv, type)=1   
                  case default             ! use -i*alm
                    a1(inv, type)=-1; a2(inv, type)= 1; swap(inv, type)=1   
                end select
            endif
        enddo
     enddo
     !$acc end kernels
#ifdef GPU
    h_swap = swap
    h_a1 = a1
    h_a2 = a2
#endif
end subroutine




subroutine get_ring_info()
#ifdef GPU
  use sht_data_module, only: nside, nring, i1_arr, i2_arr, nfft_arr, h_i1_arr, h_i2_arr, h_nfft_arr
#else
  use sht_data_module, only: nside, nring, i1_arr, i2_arr, nfft_arr
#endif
    implicit none

    integer(4)  :: i, i1, i2, nfft
    !$acc kernels
    !!$acc loop private(i1, i2, nfft)
    do i=0,nring-1
        ! compute i1, starting index of the i-th ring
        if (i.ge.0 .and. i .lt. nside) i1 = 2*i*(i+1)

        if (i.ge.nside .and. i .lt. 3*nside-1) i1 = 2*nside + 4*nside*i - 2*nside**2

        if (i.ge.3*nside-1 .and. i .lt. 4*nside-1) i1 = 12*nside**2 - 2*(4*nside-i-1)*(4*nside-i)

        ! compute nfft, size of the i-th ring
        if (i.ge.0 .and. i .lt. nside) nfft = 4*(i+1)

        if (i.ge.nside .and. i .lt. 3*nside-1) nfft = 4*nside

        if (i.ge.3*nside-1 .and. i .lt. 4*nside-1) nfft = 4*(4*nside-1-i)

        ! compute i2, end index of the i-th ring
        i2 = i1 + nfft - 1

        i1_arr(i) = i1
        i2_arr(i) = i2
        nfft_arr(i) = nfft
     enddo

     !$acc end kernels

#ifdef GPU
     h_i1_arr = i1_arr
     h_i2_arr = i2_arr
     h_nfft_arr = nfft_arr
#endif
     
end subroutine
