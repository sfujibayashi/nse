program check_HSmass
  
  use module_nuclear_data_HS
  use const, only: mpmev, mnmev, mumev
  implicit none
  
  integer :: k
  real(8) :: mass, bind, bind_check, mexc
  
  call init_nuclear_data_HS("data/dd2_frdm_comp/dd2_frdm_comp_v1.02.bin")
  
  k = find_HS_index(56,26)
  
  call get_nuclear_data_HS(k, mass, bind)
  
  bind_check = 26d0*mpmev + 30d0*mnmev - mass
  mexc = mass - 56d0*mumev
  
  write(*,'(a,es24.15)') "mass       = ", mass
  write(*,'(a,es24.15)') "bind       = ", bind
  write(*,'(a,es24.15)') "bind check = ", bind_check
  write(*,'(a,es24.15)') "mexc       = ", mexc
  write(*,'(a,es24.15)') "difference = ", bind-bind_check

end program check_HSmass
