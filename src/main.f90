program main
  use module_ptf_reaclib
  use module_nse
  use const, only : mumev
  implicit none

  real(8) :: rho, temp, ye
  integer,parameter :: itrlim = 1000
  real(8),parameter :: tol = 1d-10
!!! NSE
  integer :: n_spec
  real(8),allocatable :: xnse(:)
  integer :: k

  real(8) :: xn_history(itrlim), xp_history(itrlim)
  integer :: itr_out
  logical :: nsefail, use_tnaguess
  real(8) :: err_out

  block
    character(256) :: fn

    fn = "/Users/sfujibayashi/tmp/wanajo/winvn_v2.0.dat"
    !fn = "reduced"; nct_in = 2322; nz_in  = 55; na_in  = 90+55
    
    !call init_ptf_reaclib("/Users/fujibayashishou/Desktop/Wanajo/winvn_v2.0.dat",nct_in,nz_in,na_in)
    call init_ptf_reaclib(fn)!,nct_in,nz_in,na_in)

    !call nse_init(nct_in)
    call nse_init_reaclib(n_spec)
    ! write(6,*) "",n_spec
    ! stop
  end block

  ! call nse_init(n_spec)

  allocate(xnse(n_spec))

  use_tnaguess = .false.
  rho=1d8
  temp=6.5d9
  ye=0.30d0
  call    calc_nse(rho,temp,ye,itrlim,tol,xnse,nsefail,use_tnaguess,itr_out = itr_out, err_out = err_out)
  write(6,*) nsefail,itr_out,err_out
  if(nsefail) then
     use_tnaguess = .true.
     call calc_nse(rho,temp,ye,itrlim,tol,xnse,nsefail,use_tnaguess,itr_out = itr_out, err_out = err_out)
     write(6,*) nsefail,itr_out,err_out
  endif
  call output_composition(xnse,temp,rho,ye)
  ! call test_converge(rho,temp,ye,use_tnaguess)


  ! rho=1.6605E+03
  ! temp=1.1605E+09*2d0
  ! ye=0.30d0
  
  ! call calc_nse(rho,temp,ye,itrlim,tol,xnse,nsefail)
  block
    real(8) :: mexc_ave, z_heavy, a_heavy, y_heavy, ytot, xsum, yesum
    call statistic(xnse, mexc_ave, z_heavy, a_heavy, y_heavy, ytot, xsum, yesum)
    write(6,'("rho, temp, ye = "99es12.4)') rho,temp,ye
    write(6,'("Sanity check. sum of X_i (should be 1) and Z_i Y_i (should be Ye)  = "99es12.4)') xsum, yesum
    write(6,'("mass excess per baryon (MeV) = "99es12.4)') mexc_ave
    write(6,'(" Y = 1/<A> for all nuclei    = "99es12.4)') ytot
    write(6,'("<A>, <Z>, Y for heavy        = "99es12.4)') z_heavy, a_heavy, y_heavy
  end block

  ! call two_nuclei_approx(ye,xnse)
  ! block
  !   real(8) :: mexc_ave, z_heavy, a_heavy, y_heavy, ytot, xsum, yesum
  !   call statistic(xnse, mexc_ave, z_heavy, a_heavy, y_heavy, ytot, xsum, yesum)
  !   write(6,'(99es12.4)') rho,temp,ye, mexc_ave, z_heavy, a_heavy, y_heavy, ytot, xsum, yesum
  ! end block
  ! !call output_composition(xnse,temp,rho,ye)
  ! call test_converge(rho,temp,ye)
  ! stop

  stop

  block
    integer :: nrho=10, ntemp=10, nye=9
    integer :: irho, itemp, iye
    real(8) :: rho_max=1.66d14, rho_min=1.66d1, temp_max=10d0*1.16d10, temp_min=0.1d0*1.16d10, ye_max=0.9d0, ye_min=0.10d0
    real(8) :: mexc_ave, z_heavy, a_heavy, y_heavy, ytot, xsum, yesum
    
    integer :: k_n, k_p, k_4he
    
    k_n   = jnuc_reaclib(1,0)
    k_p   = jnuc_reaclib(1,1)
    k_4he = jnuc_reaclib(4,2)
    
    do iye=1,nye
       do itemp=1,ntemp
          write(99,*)
          do irho=1,nrho
             ye=0.5d0
             !ye = ye_min + (ye_max-ye_min)*dble(iye-1)/dble(nye-1)
             rho = 10d0**(log10(rho_min) + (log10(rho_max)-log10(rho_min))*dble(irho-1)/dble(nrho-1))
             temp = 10d0**(log10(temp_min) + (log10(temp_max)-log10(temp_min))*dble(itemp-1)/dble(ntemp-1))
             write(6,*) irho,itemp,rho,temp
             call calc_nse(rho,temp,ye,itrlim,tol,xnse,nsefail,use_tnaguess)
             call statistic(xnse, mexc_ave, z_heavy, a_heavy, y_heavy, ytot, xsum, yesum)
             write(99,'(99es12.4)') rho,temp,ye, mexc_ave/mumev, z_heavy, a_heavy, y_heavy, ytot, xnse(k_n),xnse(k_p), xnse(k_4he)
             if(nsefail)then
                write(6,'("fail",3i5,99es12.4)') irho,itemp,iye,rho,temp,ye
                stop
             endif
          enddo
       enddo
    enddo
  end block
  
end program main
