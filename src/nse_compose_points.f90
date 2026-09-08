program nse_single
  use module_nse
  use module_ptf_reaclib
  use module_ptf_rauscher
  implicit none
  
  integer,parameter :: itrlim = 300
  real(8),parameter :: tol = 1d-10

  integer :: n_spec
  real(8),allocatable :: xnse(:)

  logical :: nsefail, use_TNAguess
  
  real(8) :: mexc_ave, z_heavy, a_heavy, y_heavy, ytot, xsum, yesum, x_heavy

  integer,parameter :: n_rank = 10
  integer :: index_r(n_rank)
  integer :: i
  
  character(256) :: fn_out, fn_points

  type(stat_t) :: stat

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
    read(10,*);read(10,'(a)') fn_points
    read(10,*);read(10,'(a)') fn_out
    close(10)
    
    call init_ptf_reaclib(fn_winv)
    call init_ptf_rauscher(fn_raucher)
    call nse_init_reaclib(n_spec, use_rauscher_ptf)

  end block

  allocate(xnse(n_spec))
  
  block
    character(2048) :: line
    integer :: unit, ios
    integer :: inb, iyq, it
    integer :: npoint
    
    real(8) :: nb, temp_mev, yq
    real(8) :: rho, temp, ye
    real(8) :: a_n, z_n, y_n, abar
    real(8) :: yn, yp, yh2, yh3, yhe3, yhe4
    real(8) :: q7
    
    ! read CompOSE h5 file
    open(newunit=unit, file=trim(fn_points), status="old", action="read")

    npoint = 0

    do

       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) exit

       line = adjustl(line)

       if (len_trim(line) == 0) cycle
       if (line(1:1) == "#") cycle

       read(line,*,iostat=ios) &
            inb, iyq, it, &
            nb, temp_mev, yq, rho, temp, &
            a_n, z_n, y_n, abar, &
            yn, yp, yh2, yh3, yhe3, yhe4, q7

       if (ios /= 0) then
          write(*,*) "ERROR reading line:"
          write(*,*) trim(line)
          stop
       endif

       npoint = npoint + 1

       if (npoint == 1) then
          write(*,*) "first point:"
          write(*,'(a,3i6)')       "index       = ", inb, iyq, it
          write(*,'(a,3es16.7)')   "nb,T,Yq     = ", nb, temp_mev, yq
          write(*,'(a,2es16.7)')   "rho,T[K]    = ", rho, temp
          write(*,'(a,4es16.7)')   "A_N,Z_N,Y_N,Abar = ", &
               a_n, z_n, y_n, abar
          write(*,'(a,6es16.7)')   "Yn,Yp,Yd,Yt,YHe3,YHe4 = ", &
               yn, yp, yh2, yh3, yhe3, yhe4
          write(*,'(a,es16.7)')    "Q7          = ", q7
       endif


       if (npoint == 1) then
          ye = yq
          call calc_nse(rho,temp,ye,itrlim,tol,xnse,nsefail,use_TNAguess)
          call statistic_compose(xnse, stat)
          write(*,'(a,3es16.7)')   "rho,T[K],Ye  = ", rho, temp, ye
          write(*,'(a,4es16.7)')   "A_N,Z_N,Y_N,Abar = ", &
               stat%a_n, stat%z_n, stat%y_n, stat%abar
          write(*,'(a,6es16.7)')   "Yn,Yp,Yd,Yt,YHe3,YHe4 = ", &
               stat%yn, stat%yp, stat%yh2, stat%yh3, stat%yhe3, stat%yhe4
          stop
       endif
    enddo

    close(unit)

    write(*,*) "number of points =", npoint

  end block

  ! call calc_nse(rho,temp,ye,itrlim,tol,xnse,nsefail,use_TNAguess)
  
  ! call statistic(xnse, mexc_ave, z_heavy, a_heavy, y_heavy, ytot, xsum, yesum)
  ! x_heavy = a_heavy*y_heavy
  ! call index_rank(n_rank, n_spec, xnse, index_r, 1d0)
  
  ! write(6,'(a, 11es16.7e3)') "rho, T, Ye = ", rho,temp,ye
  ! write(6,'(10a16,10es16.7e3)') (name_reaclib(ireaclib(index_r(i))),i=1,n_rank), (xnse(index_r(i)),i=1,n_rank)

  ! call output_nse_full(rho,temp,ye,xnse, fn_out)
  
end program nse_single
