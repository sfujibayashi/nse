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
  call calc_nse(net,rho,temp,ye,itrlim,tol,xnse,nsefail,use_TNAguess)
  
  call statistic(net, xnse, mexc_ave, z_heavy, a_heavy, y_heavy, ytot, xsum, yesum)
  x_heavy = a_heavy*y_heavy

  write(6,'(a, 11es16.7e3)') "rho, T, Ye = ", rho,temp,ye
  call index_rank(n_rank, net_aprox21%n_spec, xnse_aprox21, index_r, 1d0)
  write(6,'(10a16,10es16.7e3)') (net_aprox21%name_nucl(index_r(i)),i=1,n_rank), (xnse_aprox21(index_r(i)),i=1,n_rank)

  call index_rank(n_rank, net%n_spec, xnse, index_r, 1d0)
  write(6,'(10a16,10es16.7e3)') (net%name_nucl(index_r(i)),i=1,n_rank), (xnse(index_r(i)),i=1,n_rank)

  call output_nse_full(net_aprox21, rho, temp, ye, xnse_aprox21, fn_out)
  call output_nse_full(net, rho, temp, ye, xnse, "nse_full.dat")
  
end program nse_aprox21
