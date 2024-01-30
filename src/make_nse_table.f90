program make_nse_table
  use module_nse
  use module_ptf_reaclib
  
  ! use const, only : mumev
  
  implicit none
  
  character(256) :: fn_para

  real(8) :: rho, temp, ye
  integer,parameter :: itrlim = 300
  real(8),parameter :: tol = 1d-10
!!! NSE
  integer :: n_spec
  real(8),allocatable :: xnse(:)
  integer :: k

  real(8) :: xn_history(itrlim), xp_history(itrlim)
  logical :: nsefail, use_tnaguess
  
  integer :: nrho, ntemp, nye
  real(8) :: rho_max, rho_min, temp_max, temp_min, ye_max, ye_min

  integer :: irho, itemp, iye
  integer :: itr_out
  real(8) :: err_out
  real(8) :: mexc_ave, z_heavy, a_heavy, y_heavy, ytot, xsum, yesum

  integer :: k_n, k_p, k_4he

  character(256) :: fn_out, fn_ptf

  call getarg(1, fn_para)
  
  open(10,file=fn_para,status="old",action="read")
  read(10,*);read(10,'(a)') fn_ptf
  read(10,*);read(10,'(a)') fn_out
  read(10,*);read(10,*) nrho, rho_min, rho_max
  read(10,*);read(10,*) ntemp, temp_min, temp_max
  read(10,*);read(10,*) nye, ye_min, ye_max
  close(10)

  call init_ptf_reaclib(fn_ptf)
  call nse_init_reaclib(n_spec)

  allocate(xnse(n_spec))

  k_n   = jnuc_reaclib(1,0)
  k_p   = jnuc_reaclib(1,1)
  k_4he = jnuc_reaclib(4,2)

  open(11,file=fn_out,status="replace",action="write")
  write(11,'("#",99a16)') "rho(g/cm^3)", "temp(K)", "Ye", "mexc_av(MeV)", "<Z>h", "<A>h", "Yh", "Ytot", "Xn", "Xp", "X(4He)"
  do iye=1,nye

     if(nye>1)then
        ye = ye_min + (ye_max-ye_min)*dble(iye-1)/dble(nye-1)
     else
        ye = ye_min
     endif
     
     do itemp=1,ntemp

        if(ntemp>1)then
           temp = 10d0**(log10(temp_min) + (log10(temp_max)-log10(temp_min))*dble(itemp-1)/dble(ntemp-1))
        else
           temp = temp_min
        endif

        do irho=1,nrho

           if(ntemp>1)then
              rho = 10d0**(log10(rho_min) + (log10(rho_max)-log10(rho_min))*dble(irho-1)/dble(nrho-1))
           else
              rho = rho_min
           endif
           

           use_tnaguess = .false.
           call calc_nse(rho,temp,ye,itrlim,tol,xnse,nsefail,use_tnaguess,itr_out = itr_out, err_out = err_out)
           ! write(6,*) nsefail,itr_out,err_out
           if(nsefail) then
              use_tnaguess = .true.
              call calc_nse(rho,temp,ye,itrlim,tol,xnse,nsefail,use_tnaguess,itr_out = itr_out, err_out = err_out)
              ! write(6,*) nsefail,itr_out,err_out
           endif
           write(6,*) irho,itemp,iye,rho,temp,ye,itr_out,use_tnaguess

           if(nsefail) then
              ! call calc_nse(rho,temp,ye,itrlim,tol,xnse,nsefail,use_tnaguess,itr_out = itr_out, err_out = err_out)
              call two_nuclei_approx(ye,xnse)
           endif

           call statistic(xnse, mexc_ave, z_heavy, a_heavy, y_heavy, ytot, xsum, yesum)
           write(11,'(" ",99es16.7e3)') rho,temp,ye, mexc_ave, z_heavy, a_heavy, y_heavy, ytot, xnse(k_n),xnse(k_p), xnse(k_4he)

           if(nsefail)then
              write(6,'("failed",3i5,99es12.4)') irho,itemp,iye,rho,temp,ye, err_out
              !call test_converge(rho,temp,ye,use_tnaguess)
              !call output_composition(xnse,temp,rho,ye)
              !stop
           endif


        enddo
     enddo
  enddo
  close(11)

end program make_nse_table
