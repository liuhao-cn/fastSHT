! Module of the shared data in SHT. Contain only data, no subroutine. 
!
! The initialization and freeing are controlled by the wrapper.
!
! NOTE: due to the mechanism/realization of "module variables", it is STRONGLY
! recommended and MUCH safer to first allocate the module variables using
! "sht_data_alloc_idl", which also explicitly set all variables to zero; and
! then pass external values to the module variable by "sht_set_data_idl".
!
module sht_data_module
#ifdef GPU
  use cublas_v2
  use cudafor
#endif
implicit none  
    ! number-like constants for plm-recursion
    ! diag_fac: factor for the diagonal plm 
    ! i0: start position for the ring optimization
    !real(8), save, public                                   :: mm_time = 0, fft_time = 0, plm_time = 0, tot_time = 0
    real(8),    save, public                                :: pa, pi=3.1415926535897932384626433833
    integer(4), save, public                                :: i0, i1, p1, p2, n, rlen, nhr

    ! basic SHT parameters
    integer(4), save, public                                :: nside, npix, lmax, nring, nsim, is_pol

    ! array-like constants for plm-recursion
#ifdef GPU
    real(8),    save, public, dimension(:),   allocatable, device   :: cc, ss, c1, c2, c3, c4, a, b, b_w, c_w
    real(8),    save, public, dimension(:),   allocatable, device  :: two_on_s2, normal_l, normal_m, lam_fact
    real(8),    save, public, dimension(:),   allocatable, device   :: theta, phi0, diag_fac

    ! plm0, plm1 and plm2 are slices of plm at one m, and will be used for matrix multiplication
    real(8),    save, public, dimension(:,:), allocatable, device   :: plm0, plm1, plm2

    ! the ring information
    integer(4), save, public, dimension(:),   allocatable, device   :: i1_arr, i2_arr, nfft_arr
    integer(4), save, public, dimension(:),   allocatable   :: h_i1_arr, h_i2_arr, h_nfft_arr
    
    ! factors for dgemm
    integer(4), save, public, dimension(:,:), allocatable, device   :: swap
    integer(4), save, public, dimension(:,:), allocatable   :: h_swap
    real(8),    save, public, dimension(:,:), allocatable, device   :: a1, a2
    real(8),    save, public, dimension(:,:), allocatable   :: h_a1, h_a2

    ! FFT mapping tables
    integer(4), save, public, dimension(:,:), allocatable, device   :: tab, tabI, tabJ

    ! FFT mapping factors
    real(8),    save, public, dimension(:,:), allocatable, device   :: fac1, fac2, fac3, fac4
    real(8),    save, public, dimension(:,:), allocatable, device   :: fac1I, fac2I, fac3I, fac4I

    
    ! tables for ring optimizations and safe positions
    integer(4), save, public, dimension(:),   allocatable,device   :: ring_tab

    real(8),    save, public, dimension(:,:), allocatable,device   :: val1, val2
    integer(4), save, public, dimension(:,:), allocatable,device   :: pos

    ! slices of FFT-matrix for matrix multiplication
    real(8),    save, public, dimension(:,:), allocatable, device   :: mat_cache1, mat_cache2

    ! buffs for the map-FFT
    real(8),    save, public, dimension(:,:), allocatable,pinned   :: buff1, buff2, buff3, buff4
    real(8),    save, public, dimension(:,:), allocatable, device   :: d_buff1, d_buff2, d_buff3, d_buff4
    real(8),    save, public, dimension(:,:,:), allocatable, device   :: d_alm, d_alm2
    type(cublasHandle) :: handle
    integer(cuda_stream_kind), dimension(:), allocatable :: cu_streams
    !$acc declare copyin(lmax,nring,nside,npix,pa,pi,rlen,n,i0,i1,p1,p2,nhr,nsim)

#else
    real(8),    save, public, dimension(:),   allocatable   :: cc, ss, c1, c2, c3, c4, a, b, b_w, c_w
    real(8),    save, public, dimension(:),   allocatable  :: two_on_s2, normal_l, normal_m, lam_fact
    real(8),    save, public, dimension(:),   allocatable   :: theta, phi0, diag_fac

    ! plm0, plm1 and plm2 are slices of plm at one m, and will be used for matrix multiplication
    real(8),    save, public, dimension(:,:), allocatable   :: plm0, plm1, plm2

    ! the ring information
    integer(4), save, public, dimension(:),   allocatable   :: i1_arr, i2_arr, nfft_arr
    
    ! factors for dgemm
    integer(4), save, public, dimension(:,:), allocatable   :: swap
    real(8),    save, public, dimension(:,:), allocatable   :: a1, a2

    ! FFT mapping tables
    integer(4), save, public, dimension(:,:), allocatable   :: tab, tabI, tabJ

    ! FFT mapping factors
    real(8),    save, public, dimension(:,:), allocatable   :: fac1, fac2, fac3, fac4
    real(8),    save, public, dimension(:,:), allocatable   :: fac1I, fac2I, fac3I, fac4I

    
    ! tables for ring optimizations and safe positions
    integer(4), save, public, dimension(:),   allocatable   :: ring_tab
    real(8),    save, public, dimension(:,:), allocatable   :: val1, val2
    integer(4), save, public, dimension(:,:), allocatable   :: pos

    ! slices of FFT-matrix for matrix multiplication
    real(8),    save, public, dimension(:,:), allocatable   :: mat_cache1, mat_cache2

    ! buffs for the map-FFT
    real(8),    save, public, dimension(:,:), allocatable   :: buff1, buff2, buff3, buff4
    
#endif
    ! flags for set/accumulate of the map-FFT
    integer(4), save, public, dimension(:),   allocatable   :: acc_flag1, acc_flag2

    
end module sht_data_module








