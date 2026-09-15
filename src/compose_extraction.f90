program extraction
  use module_nse
  use module_nuclear_data_winvne
  use module_ptf_rauscher
  use module_eos_helmholtz
  use module_stat_weight_policy, only: &
       stat_weight_policy_t, &
       STAT_WEIGHT_NONE, STAT_WEIGHT_WINVNE, STAT_WEIGHT_RAUSCHER, STAT_WEIGHT_HS
  use module_nse_species_policy
  
  implicit none
  
  integer,parameter :: itrlim = 300
  real(8),parameter :: tol = 1d-10

  real(8),allocatable :: xnse(:), xnse_aprox21(:)

  logical :: nsefail, use_TNAguess
  
  character(256) :: fn_out, fn_points

  type(stat_t) :: stat, stat_aprox21
  integer :: iyq_target, it_target

  real(8) :: eps, pres, cs2, entr

  type(nse_network_t) :: net, net_aprox21

  type(stat_weight_policy_t) :: stat_weight_policy
  type(nse_species_policy_t) :: species_policy

  iyq_target = 19
  it_target = 15
  
  use_TNAguess = .false.

  block
    character(256) :: fn_para,fn_helm
    character(256) :: fn_winv, fn_raucher
    logical :: use_rauscher_ptf
    
    call getarg(1, fn_para)
  
    open(10,file=fn_para,status="old",action="read")
    read(10,*);read(10,'(a)') fn_winv
    read(10,*);read(10,'(L)') use_rauscher_ptf
    read(10,*);read(10,'(a)') fn_raucher
    read(10,*);read(10,'(a)') fn_points
    read(10,*);read(10,'(a)') fn_helm
    read(10,*);read(10,'(a)') fn_out
    close(10)

    if (use_rauscher_ptf) then

       stat_weight_policy%primary  = STAT_WEIGHT_RAUSCHER
       stat_weight_policy%fallback = STAT_WEIGHT_WINVNE

    else

       stat_weight_policy%primary  = STAT_WEIGHT_WINVNE
       stat_weight_policy%fallback = STAT_WEIGHT_NONE

    endif

    species_policy%mode = NSE_SPECIES_ALL_WINVNE
    stat_weight_policy%primary  = STAT_WEIGHT_HS
    stat_weight_policy%fallback = STAT_WEIGHT_NONE
    
    
    call init_winvne(fn_winv)
    call init_ptf_rauscher(fn_raucher)
    call nse_init_winvne(net, stat_weight_policy, species_policy)
    call nse_init_aprox21(net_aprox21, stat_weight_policy)

    call init_eos(fn_helm)
  end block

  allocate(xnse(net%n_spec))
  allocate(xnse_aprox21(net_aprox21%n_spec))
  
  block
    use const, only: mu, mev2erg, mnmev, mumev
    
    character(2048) :: line
    integer :: unit, ios
    integer :: inb, iyq, it
    integer :: npoint
    
    real(8) :: nb, temp_mev, yq
    real(8) :: rho, temp, ye
    real(8) :: a_n, z_n, y_n, abar
    real(8) :: yn, yp, yh2, yh3, yhe3, yhe4
    real(8) :: q7, E_Comp_MeV

    integer :: unit_com, unit_nse, unit_a21

    real(8) :: eps_helm, pres_helm, cs2_helm, entr_helm, E_helm_MeV
    real(8) :: mres_Comp, Fcoul_Comp

    open(newunit=unit_com, file="compose.dat", status="replace", action="write")
    open(newunit=unit_nse, file="nse.dat", status="replace", action="write")
    open(newunit=unit_a21, file="nse_aprox21.dat", status="replace", action="write")
    
    write(unit_com,'("#",99a20)') "rho", "temp", "ye",  "A_N", "Z_N", "Y_N", "Abar", "Yn", "Yp", "Yh2", "Yh3", "Yhe3", "Yhe4", "E/b(MeV)", "mexc/b(with helm)", "Ecoul/b"
    write(unit_nse,'("#",99a20)') "rho", "temp", "ye",  "A_N", "Z_N", "Y_N", "Abar", "Yn", "Yp", "Yh2", "Yh3", "Yhe3", "Yhe4", "E/b(MeV)", "mexc/b", "Ecoul/b", "s/k"
    write(unit_a21,'("#",99a20)') "rho", "temp", "ye",  "A_N", "Z_N", "Y_N", "Abar", "Yn", "Yp", "Yh2", "Yh3", "Yhe3", "Yhe4", "E/b(MeV)", "mexc/b", "Ecoul/b", "s/k"
    
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

       E_Comp_MeV = mnmev*(1.d0 + q7) - mumev

       if (ios /= 0) then
          write(*,*) "ERROR reading line:"
          write(*,*) trim(line)
          stop
       endif

       ! if (mod(iyq-1,4)>0) cycle
       !if (iyq/=1.and.iyq/=9.and.iyq/=19.and.iyq/=29.and.iyq/=39.and.iyq/=49.and.iyq/=59) cycle
       !if (mod(inb-1,10)>0) cycle
       !if (mod(it-1,10)>0) cycle
       if (iyq /= iyq_target) cycle
       if (it  /= it_target)  cycle
       !if (rho > 1d14)cycle
       !if(temp < 4d9 .or. 1d10 < temp)cycle
       !if(1d12 < rho)cycle

       npoint = npoint + 1
       
       ye = yq

       call eos_all(rho, temp, ye, 1d0/abar, 0d0, &
            eps_helm, pres_helm, cs2_helm, entr_helm)

       E_helm_MeV = eps_helm * mu / mev2erg
       mres_Comp = E_Comp_MeV - E_helm_MeV
       Fcoul_Comp = y_n*fcoulomb_HS(rho, ye, z_n, a_n, net%n0_fm)
       
       call calc_nse(net,rho,temp,ye,itrlim,tol,xnse,nsefail,use_TNAguess)
       if(nsefail)then
          write(6,*) "fail in large NSE"
          error stop
       endif

       call calc_nse(net_aprox21,rho,temp,ye,itrlim,tol,xnse_aprox21,nsefail,use_TNAguess)
       if(nsefail)then
          call calc_nse(net_aprox21,rho,temp,ye,itrlim,tol,xnse_aprox21,nsefail,.true.)
       endif
       if(nsefail)then
          write(6,*) "fail in aprox21"
          error stop
       endif
       call statistic_compose(net, rho, xnse, stat)
       call statistic_compose(net_aprox21, rho, xnse_aprox21, stat_aprox21)
       write(unit_com,'(" ",99es20.11e3)') rho, temp, ye,  a_n, z_n, y_n, abar, yn, yp, yh2, yh3, yhe3, yhe4, E_Comp_MeV, mres_Comp, Fcoul_Comp
       
       call eos_all(rho, temp, ye, 1d0/stat%abar, stat%mexc, &
            eps, pres, cs2, entr)
            
       write(unit_nse,'(" ",99es20.11e3)') rho, temp, ye, stat%a_n, stat%z_n, stat%y_n, stat%abar, stat%yn, stat%yp, stat%yh2, stat%yh3, stat%yhe3, stat%yhe4, E_helm_MeV, stat%mexc, stat%ecoul, entr

       call eos_all(rho, temp, ye, 1d0/stat_aprox21%abar, stat_aprox21%mexc, &
            eps, pres, cs2, entr)
            
       write(unit_a21,'(" ",99es20.11e3)') rho, temp, ye, stat_aprox21%a_n, stat_aprox21%z_n, stat_aprox21%y_n, stat_aprox21%abar, stat_aprox21%yn, stat_aprox21%yp, stat_aprox21%yh2, stat_aprox21%yh3, stat_aprox21%yhe3, stat_aprox21%yhe4, E_helm_MeV, stat_aprox21%mexc, stat_aprox21%ecoul, entr
       
       write(6,*) inb,temp,rho,yq

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
  
end program extraction
