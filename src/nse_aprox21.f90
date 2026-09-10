program nse_aprox21
  use module_nse
  use module_ptf_reaclib
  use module_ptf_rauscher
  implicit none
  
  real(8) :: rho, temp, ye
  integer,parameter :: itrlim = 300
  real(8),parameter :: tol = 1d-10

  type(nse_network_t) :: net_aprox21, net
  real(8),allocatable :: xnse(:), xnse_aprox21(:)
  type(stat_t) :: stat, stat_aprox21

  logical :: nsefail, use_TNAguess
  
  real(8) :: mexc_ave, z_heavy, a_heavy, y_heavy, ytot, xsum, yesum, x_heavy

  integer,parameter :: n_rank = 10
  integer :: index_r(n_rank)
  integer :: i
  
  character(256) :: fn_out

  use_TNAguess = .false.

  block
    character(256) :: fn_para
    character(256) :: fn_winv, fn_raucher
    logical :: use_rauscher_ptf
    
    call getarg(1, fn_para)
  
    open(10,file=fn_para,status="old",action="read")
    read(10,*);read(10,'(a)') fn_winv
    read(10,*);read(10,'(L)') use_rauscher_ptf
    read(10,*);read(10,'(a)') fn_raucher
    read(10,*);read(10,'(a)') fn_out
    read(10,*);read(10,*) rho, temp, ye
    close(10)
    
    call init_ptf_reaclib(fn_winv)
    call init_ptf_rauscher(fn_raucher)
    call nse_init_reaclib(net, use_rauscher_ptf)
    call nse_init_aprox21(net_aprox21, use_rauscher_ptf)

  end block

  allocate(xnse_aprox21(net_aprox21%n_spec))
  allocate(xnse(net%n_spec))

  call calc_nse(net_aprox21,rho,temp,ye,itrlim,tol,xnse_aprox21,nsefail,use_TNAguess)
  if(nsefail)then
     write(6,*) "NSE does not converge in aprox21"
  endif
  call calc_nse(net,rho,temp,ye,itrlim,tol,xnse,nsefail,use_TNAguess)
  if(nsefail)then
     write(6,*) "NSE does not converge in large set"
  endif
  
  call output_nse_full(net_aprox21, rho, temp, ye, xnse_aprox21, fn_out)
  call output_nse_full(net, rho, temp, ye, xnse, "nse_full.dat")

  write(6,'(a, 11es16.7e3)') "rho, T, Ye = ", rho,temp,ye
  call statistic_compose(net, rho, xnse, stat)
  call statistic_compose(net_aprox21, rho, xnse_aprox21, stat_aprox21)
  ! call index_rank(n_rank, net_aprox21%n_spec, xnse_aprox21, index_r, 1d0)
  ! write(6,'(10a16,10es16.7e3)') (net_aprox21%name_nucl(index_r(i)),i=1,n_rank), (xnse_aprox21(index_r(i)),i=1,n_rank)

  ! call index_rank(n_rank, net%n_spec, xnse, index_r, 1d0)
  ! write(6,'(10a16,10es16.7e3)') (net%name_nucl(index_r(i)),i=1,n_rank), (xnse(index_r(i)),i=1,n_rank)


  write(6,'(a,99es20.11e3)') "aprox21: ",stat_aprox21%a_n, stat_aprox21%z_n, stat_aprox21%y_n, stat_aprox21%abar, stat_aprox21%yn, stat_aprox21%yp, stat_aprox21%yh2, stat_aprox21%yh3, stat_aprox21%yhe3, stat_aprox21%yhe4, stat_aprox21%mexc, stat_aprox21%ecoul
  write(6,'(a,99es20.11e3)') "full   : ",stat%a_n, stat%z_n, stat%y_n, stat%abar, stat%yn, stat%yp, stat%yh2, stat%yh3, stat%yhe3, stat%yhe4, stat%mexc, stat%ecoul

  
end program nse_aprox21
