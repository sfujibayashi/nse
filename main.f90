program main
  use module_ptf_reaclib
  use module_nse
  implicit none

  real(8) :: rho, temp, ye
  integer,parameter :: itrlim = 200
  real(8),parameter :: tol = 1d-10
!!! NSE
  integer :: n_spec
  real(8),allocatable :: xnse(:)
  integer :: k

  real(8) :: xn_history(itrlim), xp_history(itrlim)
  integer :: itr_out

  block
    character(256) :: fn
    integer :: nct_in, nz_in, na_in

    fn = "/Users/fujibayashishou/Desktop/Wanajo/winvn_v2.0.dat"; nct_in = 7854; nz_in  = 112; na_in  = 337
    !fn = "reduced"; nct_in = 2322; nz_in  = 55; na_in  = 90+55
    
    !call init_ptf_reaclib("/Users/fujibayashishou/Desktop/Wanajo/winvn_v2.0.dat",nct_in,nz_in,na_in)
    call init_ptf_reaclib(fn,nct_in,nz_in,na_in)

    !call nse_init(nct_in)
    call nse_init_reaclib(n_spec)
    write(6,*) n_spec
    !stop
  end block

  allocate(xnse(n_spec))

  rho=1d8
  temp=1d9
  ye=0.50d0
  call calc_nse(rho,temp,ye,itrlim,tol,xnse)
  call output_composition(xnse,temp,rho,ye)
  
  !call test_converge(rho,temp,ye)
  
  ! do k=1,n_spec
  !    write(99,'(2i5,99es12.4)') nnt_reaclib(k), npt_reaclib(k), max(1d-99,xnse(k))
  ! enddo
  
end program main
