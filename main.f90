program main
  !use module_nuclei_quantum
  use module_nse
  implicit none

  real(8) :: rho, temp, ye
  integer,parameter :: itrlim = 200
  real(8),parameter :: tol = 1d-15
!!! NSE
  integer :: n_spec
  real(8),allocatable :: xnse(:)
  
  call nse_init(n_spec)

  allocate(xnse(n_spec))

  rho=1d8
  temp=10d9
  ye=0.34d0
  call calc_nse(rho,temp,ye,itrlim,tol,xnse)
  
end program main
