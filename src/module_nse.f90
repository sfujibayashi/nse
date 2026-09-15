module module_nse
  use module_stat_weight_policy, only: stat_weight_policy_t, stat_weight_ref_t, resolve_stat_weight_ref, &
       STAT_WEIGHT_NONE, STAT_WEIGHT_WINVNE, STAT_WEIGHT_RAUSCHER, STAT_WEIGHT_HS, &
       valid_stat_weight_policy

  implicit none

  private
  public :: nse_init_four,nse_init_aprox21,nse_init_winvne
  public :: calc_nse,test_converge,output_composition,statistic,statistic_compose, two_nuclei_approx

  public :: output_nse_full
  public :: fcoulomb_HS
  
  
  type :: stat_t
    real(8) :: yn
    real(8) :: yp
    real(8) :: yh2
    real(8) :: yh3
    real(8) :: yhe3
    real(8) :: yhe4
    real(8) :: a_n
    real(8) :: z_n
    real(8) :: y_n
    real(8) :: abar
    real(8) :: mexc
    real(8) :: ecoul
 end type stat_t

 public :: stat_t
 
 type :: nse_network_t
    integer :: n_spec = 0
    
    real(8), allocatable :: mexc(:)
    real(8), allocatable :: a(:)
    real(8), allocatable :: z(:)
    real(8), allocatable :: n(:)
    real(8), allocatable :: zai(:)
    real(8), allocatable :: g0(:)
    
    character(5), allocatable :: name_nucl(:)
    
    logical :: use_winvne = .false.
    
    real(8) :: n0_fm = 0.1583d0

    type(stat_weight_policy_t) :: stat_weight_policy
    type(stat_weight_ref_t), allocatable :: stat_weight(:)
    
 end type nse_network_t

  public :: nse_network_t

contains
  
  subroutine nse_init_four(net)
    use const,only:mnmev,mpmev,mamev,mumev,memev
    type(nse_network_t), intent(out) :: net
    real(8),parameter :: mexc_56ni_mev = -53.907539d0

        
    net%use_winvne = .false.

    net%n_spec = 4
    allocate(net%mexc(net%n_spec), net%a(net%n_spec), net%z(net%n_spec), net%n(net%n_spec), net%g0(net%n_spec), net%zai(net%n_spec))
    ! n
    net%a(1) = 1d0; net%z(1) = 0d0; net%n(1) = 1d0; net%g0(1) = 2d0; net%mexc(1) = mnmev-net%a(1)*mumev
    ! p
    net%a(2) = 1d0; net%z(2) = 1d0; net%n(2) = 0d0; net%g0(2) = 2d0; net%mexc(2) = mpmev-net%a(2)*mumev
    ! alpha
    net%a(3) = 4d0; net%z(3) = 2d0; net%n(3) = 2d0; net%g0(3) = 1d0; net%mexc(3) = mamev-net%a(3)*mumev
    ! 56Ni
    net%a(4) =56d0; net%z(4) =28d0; net%n(4) =28d0; net%g0(4) = 1d0; net%mexc(4) = mexc_56ni_mev-net%z(4)*memev

    net%zai(:) = net%z(:)/net%a(:)
    
  end subroutine nse_init_four

  subroutine nse_init_aprox21(net, stat_weight_policy)

    use module_nuclear_data_winvne
    use module_ptf_rauscher, only: find_rauscher_index
    use const, only: memev

    type(nse_network_t), intent(out) :: net
    type(stat_weight_policy_t),intent(in) :: stat_weight_policy

    integer, parameter :: ns = 20
    integer, parameter :: aa(ns) = [ &
         1, 1, 3, 4, 12, 14, 16, 20, 24, 28, &
         32,36,40,44,48,56,52,54,56,56 ]
    integer, parameter :: zz(ns) = [ &
         0, 1, 2, 2,  6,  7,  8, 10, 12, 14, &
         16,18,20,22,24,24,26,26,26,28 ]

    integer :: i

    integer, allocatable :: iwinvne(:), irauscher(:)
    
    allocate(iwinvne(ns), irauscher(ns))
    
    net%n_spec = ns
    net%use_winvne = .true.

    if (.not. valid_stat_weight_policy(stat_weight_policy)) then
       write(*,*) "ERROR: invalid statistical-weight policy"
       error stop
    endif
    
    net%stat_weight_policy = stat_weight_policy

    allocate(net%name_nucl(ns))
    allocate(net%mexc(ns), net%a(ns), net%z(ns), net%n(ns), &
         net%g0(ns), net%zai(ns))

    do i = 1, ns

       iwinvne(i) = find_winvne_index(aa(i), zz(i))
       irauscher(i) = find_rauscher_index(aa(i), zz(i))

       if (iwinvne(i) <= 0) then
          write(*,*) "ERROR: aprox21 nucleus missing from WinVNE:", aa(i), zz(i)
          stop
       endif

       net%a(i) = ams_winvne(iwinvne(i))
       net%z(i) = dble(npt_winvne(iwinvne(i)))
       net%n(i) = dble(nnt_winvne(iwinvne(i)))

       net%mexc(i) = exc_winvne(iwinvne(i)) - net%z(i)*memev
       net%name_nucl(i) = name_winvne(iwinvne(i))

    enddo

    net%zai(:) = net%z(:)/net%a(:)

    call resolve_network_stat_weights( &
         net, iwinvne, irauscher)
    
    deallocate(iwinvne, irauscher)

  end subroutine nse_init_aprox21

  subroutine nse_init_winvne(net, stat_weight_policy)

    use module_nuclear_data_winvne
    use module_ptf_rauscher, only: find_rauscher_index
    use const, only: memev
    
    type(nse_network_t),intent(out) :: net
    type(stat_weight_policy_t), intent(in) :: stat_weight_policy

    integer :: i, k
    integer,allocatable :: jrauscher(:)

    integer, allocatable :: iwinvne(:), irauscher(:)

    if (.not. valid_stat_weight_policy(stat_weight_policy)) then
       write(*,*) "ERROR: invalid statistical-weight policy"
       error stop
    endif
    
    net%stat_weight_policy = stat_weight_policy

    allocate(jrauscher(nct_winvne))
    jrauscher(:) = 0
    
    ! ---------------------------------------------------------
    ! WinVNE index -> Rauscher index
    ! ---------------------------------------------------------

    do k=1,nct_winvne

       jrauscher(k) = find_rauscher_index(naw_winvne(k), npt_winvne(k))

    enddo

    ! ---------------------------------------------------------
    ! Keep:
    !   1. nuclei present in Rauscher
    !   2. light nuclei missing from Rauscher -> WinVNE fallback
    !
    ! Exclude:
    !   Rauscher-missing nuclei with Z >= 87
    ! ---------------------------------------------------------

    net%n_spec = count((jrauscher > 0) .or. (npt_winvne < 87))

    write(6,'(a,i6)') "Rauscher matched species = ", &
         count(jrauscher > 0)

    write(6,'(a,i6)') "WinVNE PF fallback       = ", &
         count((jrauscher == 0) .and. (npt_winvne < 87))

    write(6,'(a,i6)') "Excluded species         = ", &
         count((jrauscher == 0) .and. (npt_winvne >= 87))

    write(6,'(a,i6)') "NSE species kept         = ", net%n_spec


    allocate(iwinvne(net%n_spec), irauscher(net%n_spec))

    allocate(net%name_nucl(net%n_spec))
    allocate(net%mexc(net%n_spec), net%a(net%n_spec), net%z(net%n_spec), net%n(net%n_spec), &
         net%g0(net%n_spec), net%zai(net%n_spec))

    ! ---------------------------------------------------------
    ! Construct NSE species arrays
    ! ---------------------------------------------------------

    i = 0

    do k=1,nct_winvne

       if (jrauscher(k) > 0 .or. npt_winvne(k) < 87) then

          i = i + 1
          
          iwinvne(i)  = k
          irauscher(i) = jrauscher(k)

          net%a(i) = ams_winvne(k)
          net%z(i) = dble(npt_winvne(k))
          net%n(i) = dble(nnt_winvne(k))

          net%mexc(i) = exc_winvne(k) - net%z(i)*memev
          net%name_nucl(i) = name_winvne(k)

       endif

    enddo

    if (i /= net%n_spec) then
       write(*,*) "ERROR constructing NSE species:", i, net%n_spec
       stop
    endif

    net%zai(:) = net%z(:)/net%a(:)

    net%use_winvne = .true.

    deallocate(jrauscher)

    call resolve_network_stat_weights( &
         net, iwinvne, irauscher)
    deallocate(iwinvne, irauscher)
    
  end subroutine nse_init_winvne


  subroutine calc_coulomb(net,rho,ye,fcoul)
    use const, only : qe, mu, pi, mev2erg

    type(nse_network_t), intent(in) :: net
    real(8),intent(in) :: rho,ye
    real(8),intent(out) :: fcoul(net%n_spec)
    real(8) :: n0
    integer :: k
    real(8) :: ne, v_n, v_c, u
    
    n0 = 0.16d0*1d39
    
    fcoul(:) = 0d0
    
    ne = ye*rho/mu
    do k=1,net%n_spec
       if( net%z(k)>0.d0)then
          v_n = net%a(k)/n0
          v_c = net%z(k)/ne
          u   = v_n/v_c
          
          fcoul(k) = (3d0/5d0)*(4d0*pi/3d0)**(-1d0/3d0) * qe**2 * n0**2 * (net%z(k)/net%a(k))**2 * (v_n)**(5d0/3d0) &
               * (-3d0/2d0*u**(1d0/3d0) + 1d0/2d0*u) / mev2erg
       endif
    enddo
    
  end subroutine calc_coulomb

  subroutine calc_ptf_nse(net, t9, g)

    use module_nuclear_data_winvne, only: get_stat_weight_winvne
    use module_ptf_rauscher, only: get_stat_weight_rauscher
    use module_stat_weight_policy, only: &
         STAT_WEIGHT_NONE, STAT_WEIGHT_WINVNE, STAT_WEIGHT_RAUSCHER
    use module_stat_weight_HS, only: get_stat_weight_HS

    type(nse_network_t), intent(in) :: net
    real(8),intent(in)  :: t9
    real(8),intent(out) :: g(net%n_spec)

    integer :: i
    
    do i = 1, net%n_spec

       select case (net%stat_weight(i)%source)
          
       case (STAT_WEIGHT_WINVNE)
          
          call get_stat_weight_winvne( &
               t9, net%stat_weight(i)%index, g(i))
          
       case (STAT_WEIGHT_RAUSCHER)
          
          call get_stat_weight_rauscher( &
               t9, net%stat_weight(i)%index, g(i))

       case (STAT_WEIGHT_HS)
          
          call get_stat_weight_HS( &
               t9, net%a(i), net%z(i), net%n(i), net%mexc(i), g(i))
          
       case default
          
          write(*,*) "ERROR: unresolved statistical weight for ", &
               net%name_nucl(i)
          error stop
          
       end select

    enddo

  end subroutine calc_ptf_nse


  subroutine calc_coulomb_HS(net, rho, ye, fcoul)

    use const, only : mu

    type(nse_network_t), intent(in) :: net
    real(8),intent(in)  :: rho, ye
    real(8),intent(out) :: fcoul(net%n_spec)

    real(8) :: ne, n0
    integer :: k

    ! fm^-3 -> cm^-3
    n0 = net%n0_fm * 1d39

    ne = ye*rho/mu

    fcoul(:) = 0d0

    do k=1,net%n_spec
       fcoul(k) = fcoulomb_HS(rho, ye, net%z(k), net%a(k), net%n0_fm)
    enddo
    
  end subroutine calc_coulomb_HS
  
  function fcoulomb_HS(rho, ye, z, a, n0_fm) result(fcoul)

    use const, only : mu, pi, fine, hbar, clight, mev2erg
    
    real(8),intent(in)  :: rho, ye, z, a, n0_fm
    real(8) :: fcoul
    real(8) :: ne, n0, r, x

    ! fm^-3 -> cm^-3
    n0 = n0_fm * 1d39

    ne = ye*rho/mu

    ! no Coulomb correction for free proton
    if (a <= 1d0 .or. z <= 0d0)then
       fcoul = 0d0
       return
    endif
    
    r = (3d0*a/(4d0*pi*n0))**(1d0/3d0)
    
    x = (ne/n0 * a/z)**(1d0/3d0)
    
    fcoul = -3d0/5d0 * z**2 * fine*hbar*clight/r &
         * (1.5d0*x - 0.5d0*x**3) / mev2erg
    
  end function fcoulomb_HS
  
  
  integer function find_nucleus(net, ia, iz) result(idx)
    type(nse_network_t), intent(in) :: net
    integer, intent(in) :: ia, iz
    integer :: i
    
    idx = 0
    do i = 1, net%n_spec
       if (nint(net%a(i)) == ia .and. nint(net%z(i)) == iz) then
          idx = i
          return
       endif
    enddo
  end function find_nucleus

  
  subroutine calc_nse(net, rho,temp,ye,itrlim,tol,xnse,nsefail,use_TNAguess,xn_history,xp_history,itr_out,err_out,xn_guess,xp_guess,xn_out,xp_out)
    use const,only : mu,kerg,pi,hbar,mev2erg
    use module_nuclear_data_winvne
    type(nse_network_t),intent(in) :: net
    real(8),intent(in) :: rho,temp,ye
    integer,intent(in) :: itrlim
    real(8),intent(in) :: tol
    real(8),intent(out) :: xnse(net%n_spec)
    logical,intent(out) :: nsefail
    logical,intent(in) :: use_TNAguess
    real(8),intent(out),optional :: xn_history(0:itrlim),xp_history(0:itrlim)
    integer,intent(out),optional :: itr_out
    real(8),intent(out),optional :: err_out
    real(8),intent(in),optional :: xn_guess,xp_guess
    real(8),intent(out),optional :: xn_out,xp_out
    
    real(8) :: logrho0
    real(8) :: logge(net%n_spec), logx(net%n_spec), fcoul(net%n_spec), g(net%n_spec)
    
    real(8) :: xp,xn

    integer :: itr
    !integer,parameter :: itrlim=50
    !    real(8),parameter :: tol = 1d-13
    
    real(8) :: dxdp,dxdn,dyedp,dyedn,dx,dye,det,dxn,dxp,dl,fac

    real(8) :: t9

    real(8),parameter :: n0 = 0.16d0*1d39
    
    ! log(rho0/rho)
    logrho0 = 2.5d0*log(mu) + 1.5d0*log(kerg*temp) - 1.5d0*log(2d0*pi) - 3d0*log(hbar) - log(rho)
    
    ! partition function may be calculated here
    t9 = temp/1d9
    if(net%use_winvne)then
       call calc_ptf_nse(net, t9,g)
    else
       g(:) = net%g0(:)
    endif
    
    if(net%use_winvne)then
       call calc_coulomb_HS(net, rho, ye, fcoul)
    else
       fcoul(:) = 0d0
    endif
    !
    logge(:) = log(g(:)) + 2.5d0*log(net%a(:)) + logrho0 - net%mexc(:)*mev2erg/(kerg*temp) &
         - fcoul(:)*mev2erg/(kerg*temp)
    
    if(present(xn_history)) xn_history(:) = 0d0
    if(present(xp_history)) xp_history(:) = 0d0

    ! use two-nuclei approx for initial guess
    if(use_TNAguess)then
       block
         integer :: k1,k2
         real(8) :: eta01ex, eta02ex
         real(8) :: z1,z2,a1,a2,x1,x2,g1,g2,mex1,mex2,n1,n2

         if(ye/=0.5d0)then
            call two_nuclei_approx_index(net, ye, fcoul, k1, k2)
         else
            k1 = find_nucleus(net, 1, 1)
            k2 = find_nucleus(net, 56, 26)
            if (k1 == 0 .or. k2 == 0) then
               write(*,*) "ERROR in calc_nse: required nucleus not found"
               write(*,*) "p index    =", k1
               write(*,*) "Fe56 index =", k2
               error stop
            endif
         endif

         z1 = net%z(k1)
         z2 = net%z(k2)
         a1 = net%a(k1)
         a2 = net%a(k2)
         n1 = a1-z1
         n2 = a2-z2
         x2 = (z1/a1 - ye)/(z1/a1 - z2/a2)
         x1 = (1d0 - x2)
         
         g1 = g(k1)
         g2 = g(k2)
         mex1 = net%mexc(k1)*mev2erg
         mex2 = net%mexc(k2)*mev2erg
         
         ! (mu_1 - m_1 c^2 + mexc_1*c^2)/kT
         eta01ex = -logrho0 + log(x1) - log(g1) - 2.5d0*log(a1) + (net%mexc(k1) + fcoul(k1))*mev2erg/(kerg*temp)
         ! (mu_2 - m_2 c^2 + mexc_2*c^2)/kT
         eta02ex = -logrho0 + log(x2) - log(g2) - 2.5d0*log(a2) + (net%mexc(k2) + fcoul(k2))*mev2erg/(kerg*temp)
         
         xn = (z2*eta01ex - z1*eta02ex)/(n1*z2-n2*z1)
         xp = (n2*eta01ex - n1*eta02ex)/(n2*z1-n1*z2)

         ! write(6,*) z1,a1,z2,a2,x1,x2, xn, xp
         ! write(6,*) (z2*eta01ex - z1*eta02ex)/(n1*z2-n2*z1), (n2*eta01ex - n1*eta02ex)/(n2*z1-n1*z2)

       end block

    else

       if(present(xn_guess).and.present(xp_guess))then
          xp = xp_guess
          xn = xn_guess
       else
          xn = -2d0*log(10d0) - logge(1)! + 100d0*ye/(temp/1.16d9)*0d0
          xp = -2d0*log(10d0) - logge(2)! - 100d0*ye/(temp/1.16d9)*0d0
          do itr=1,10000
             call step(net,xn,xp,ye,logge,dx,dye,dxn,dxp,det,dxdn,dxdp,dyedn,dyedp)
             logx(:) = logge(:) + net%z(:)*xp + net%n(:)*xn
             !write(6,'(99es12.4)') xn,xp,det,maxval(logx(:)),minval(logx(:))
             if( abs(det)>0d0 .and. dx<0d0 .and.dye<0d0)then
                exit
             else
                !xn = xn - 1d0
                !xp = xp - 1d0
                xn = xn-log(10d0)
                xp = xp-log(10d0)
             endif
             
          enddo
       endif

    endif

    !write(6,'(99es12.4)') xn,xp
    !stop
    ! call step(xn,xp,ye,logge,dx,dye,dxn,dxp,det,dxdn,dxdp,dyedn,dyedp)
    ! !write(6,'(99es12.4)') xn,xp,ye,logge,dx,dye,dxn,dxp,det,dxdn,dxdp,dyedn,dyedp

    ! if(abs(det)/min(abs(dxdn),abs(dxdp),abs(dyedn),abs(dyedp))<1d-15)then
    !    stop "not converged 0"

    !    ! xnse(3) = min(ye,1d0-ye)*2d0
    !    ! xnse(1) = max(1d-99, 1d0-ye - 0.5d0*xnse(3))
    !    ! xnse(2) = max(1d-99, ye     - 0.5d0*xnse(3))
    !    ! if( abs(ye-0.5d0)<1d-16 )then
    !    !    xnse(1) = 1d-99
    !    !    xnse(2) = 1d-99
    !    !    xnse(3) = 1d0
    !    ! endif
    !    ! if(present(itr_out)) itr_out = 0
    !    ! if(present(xn_history)) xn_history(0) = 0d0
    !    ! if(present(xp_history)) xp_history(0) = 0d0
    !    ! return
    ! endif

    nsefail = .false.
    
    if(present(xn_history)) xn_history(0) = xn
    if(present(xp_history)) xp_history(0) = xp
    if(present(itr_out))itr_out = 0
    
    do itr=1,itrlim

       if(present(xn_history))xn_history(itr) = xn
       if(present(xp_history))xp_history(itr) = xp
       if(present(itr_out))itr_out = itr

       call step(net,xn,xp,ye,logge,dx,dye,dxn,dxp,det,dxdn,dxdp,dyedn,dyedp)
       
       if(present(err_out))err_out = max(abs(dx), abs(dye))

       
       !if( (abs(dx)<tol .and. abs(dye) < tol) .or. (abs(dxn/xn)<tol .and. abs(dxp/xp)<tol) ) exit
       if( abs(dx) < tol .and. abs(dye) < tol ) exit
       
       ! if( abs(det)/min(abs(dxdn),abs(dxdp),abs(dyedn),abs(dyedp))<1d-15 .or. logge(1) + xn < -3d2 .or. logge(2) + xp < -3d2 )then
       !    ! xnse(3) = min(ye,1d0-ye)*2d0
       !    ! xnse(1) = max(1d-99, 1d0-ye - 0.5d0*xnse(3))
       !    ! xnse(2) = max(1d-99, ye     - 0.5d0*xnse(3))
       !    ! if( abs(ye-0.5d0)<1d-16 )then
       !    !    xnse(1) = 1d-99
       !    !    xnse(2) = 1d-99
       !    !    xnse(3) = 1d0
       !    ! endif
       !    write(6,*) "not converged"
       !    nsefail = .true.
       !    return
       ! endif
       
       dl = sqrt(dxn*dxn+dxp*dxp)

       ! write(6,'(i5,99es15.7)') itr,xn,xp,dx,dye,dxn,dxp,det,dxdn,dxdp,dyedn,dyedp,dl
       fac = 1d0
       if(dl>0.5d0*log(10d0))fac = 0.5d0*log(10d0)/dl
       !if(max(abs(dxn/xn),abs(dxp/xp)) > 0.5d0) fac = 0.5d0/max(abs(dxn/xn),abs(dxp/xp))
       !write(6,*) xn,xp
       
       xn = xn + dxn*fac
       xp = xp + dxp*fac

       !itr_out = itr
       if(itr==itrlim)then
          nsefail = .true.
       endif
    enddo

    logx(:) = logge(:) + net%z(:)*xp + net%n(:)*xn

    block
      real(8) :: logx_max, u(net%n_spec)
      logx_max = maxval(logx(:))
      u(:) = exp(logx(:) - logx_max)
      xnse(:) = u(:) / sum(u(:))
    end block
    if(present(xn_out)) xn_out = xn
    if(present(xp_out)) xp_out = xp

  end subroutine calc_nse

  subroutine two_nuclei_approx(net,rho,ye,xnse)
    type(nse_network_t),intent(in) :: net
    real(8),intent(in) :: rho,ye
    real(8),intent(out) :: xnse(net%n_spec)
    
    real(8) :: y1,y2, z1,z2,a1,a2
    integer :: k1,k2

    real(8) :: fcoul(net%n_spec)
    
    if(net%use_winvne)then
       call calc_coulomb_HS(net, rho, ye, fcoul)
    else
       fcoul(:) = 0d0
    endif

    call two_nuclei_approx_index(net, ye, fcoul, k1, k2)

    z1 = net%z(k1)
    z2 = net%z(k2)
    a1 = net%a(k1)
    a2 = net%a(k2)
    y2 = (z1/a1 - ye)/(z1/a1 - z2/a2)/a2
    y1 = (1d0 - a2*y2)/a1
    
    xnse(:) = 1d-300
    xnse(k1) = a1*y1
    xnse(k2) = a2*y2
    
  end subroutine two_nuclei_approx


  subroutine two_nuclei_approx_index(net,ye,fcoul,k1_min,k2_min)
    type(nse_network_t),intent(in) :: net
    real(8),intent(in) :: ye, fcoul(net%n_spec)
    integer,intent(out) :: k1_min,k2_min
    real(8) :: y1,y2, z1,z2,a1,a2, mexc1, mexc2, f, f_min

    integer :: k1, k2

    k1_min = 0
    k2_min = 0
    f_min = 1d99
    do k1=1,net%n_spec
       do k2=1,k1-1
          z1 = net%z(k1)
          z2 = net%z(k2)
          a1 = net%a(k1)
          a2 = net%a(k2)
          if( (z1/a1 - ye)*(z1/a1 - z2/a2) > 0d0 )then
             y2 = (z1/a1 - ye)/(z1/a1 - z2/a2)/a2
             y1 = (1d0 - a2*y2)/a1
             
             mexc1 = net%mexc(k1)
             mexc2 = net%mexc(k2)
             
             f = (mexc1+fcoul(k1))*y1 + (mexc2+fcoul(k2))*y2

             !write(6,*) k1,k2,y1,y2,f
             if(y1>0d0 .and. f < f_min)then
                f_min = f
                k1_min = k1
                k2_min = k2
                
                !write(6,'(2i5,99es12.4)') k1_min,k2_min,f_min,mexc1,mexc2, y1,y2
                !stop
             endif
          endif
       enddo
    enddo
    
  end subroutine two_nuclei_approx_index

  subroutine test_converge(net,rho,temp,ye,use_TNAguess)
    use const,only : mu,kerg,pi,hbar,mev2erg
    use module_nuclear_data_winvne
    type(nse_network_t),intent(in) :: net
    real(8),intent(in) :: rho,temp,ye
    logical,intent(in) :: use_TNAguess
    real(8) :: logrho0
    real(8) :: logge(net%n_spec), fcoul(net%n_spec)
    
    real(8) :: xp,xn

    real(8) :: dxdp,dxdn,dyedp,dyedn,dx,dye,det,dxn,dxp

    integer :: in,ip,nn,np
    real(8) :: xn_min,xn_max,xp_min,xp_max

    integer :: itr,itr_out
    real(8) :: tol = 1d-12
    integer,parameter :: itrlim=1000
    real(8) :: xn_history(0:itrlim),xp_history(0:itrlim),xnse(net%n_spec)

    ! real(8) :: xm,dxm,xm_min,xm_max
    logical :: nsefail
    real(8) :: t9
    real(8) :: g(net%n_spec)

    call calc_nse(net,rho,temp,ye,itrlim,tol,xnse,nsefail,use_TNAguess,xn_history,xp_history,itr_out)
    !write(6,*) xnse(:)
    write(6,*) itr_out

    open(98,file="convergence_history.dat",status="replace",action="write")
    do itr = 0,itr_out
       write(98,*) itr,xn_history(itr),xp_history(itr)
    enddo
    close(98)
    
    logrho0 = 2.5d0*log(mu) + 1.5d0*log(kerg*temp) - 1.5d0*log(2d0*pi) - 3d0*log(hbar) - log(rho)

    ! partition function may be calculated here
    t9 = temp/1d9
    if(net%use_winvne)call calc_ptf_nse(net,t9,g)

    if(net%use_winvne)then
       call calc_coulomb(net,rho,ye,fcoul)
    else
       fcoul(:) = 0d0
    endif
    !
    logge(:) = log(g(:)) + 2.5d0*log(net%a(:)) + logrho0 - net%mexc(:)*mev2erg/(kerg*temp) &
         - fcoul(:)*mev2erg/(kerg*temp)

    nn=100
    np=100
    xn_min = -10d0
    xn_max =  10d0
    xp_min = -10d0
    xp_max =  10d0

    if(xn_min>xn_history(itr_out)) xn_min = xn_history(itr_out)-5d0
    if(xp_min>xp_history(itr_out)) xp_min = xp_history(itr_out)-5d0
    if(xn_max<xn_history(itr_out)) xn_max = xn_history(itr_out)+5d0
    if(xp_max<xp_history(itr_out)) xp_max = xp_history(itr_out)+5d0

    dxn = maxval(xn_history(1:itr_out)) - minval(xn_history(1:itr_out))
    dxp = maxval(xp_history(1:itr_out)) - minval(xp_history(1:itr_out))

    xn_min = minval(xn_history(1:itr_out))-dxn*0.1d0
    xn_max = maxval(xn_history(1:itr_out))+dxn*0.1d0
    xp_min = minval(xp_history(1:itr_out))-dxp*0.1d0
    xp_max = maxval(xp_history(1:itr_out))+dxp*0.1d0
    !xm_min = -20d0
    !xm_max = 20d0

    
    open(99,file="convergence_map.dat",status="replace",action="write")
    write(99,'("#",99es12.4)') rho,temp,ye
    write(99,'("#",99es12.3e3)') xnse(1:3)
    
    do ip=1,np
       write(99,*)
       do in=1,nn

          xn = xn_min + (xn_max-xn_min)*dble(in-1)/dble(nn-1)
          xp = xp_min + (xp_max-xp_min)*dble(ip-1)/dble(np-1)

          !write(6,*) in,ip,xn,xp

          !xm = xm_min + (xm_max-xm_min)*dble(in-1)/dble(nn-1)

          !call step2(xp,xm,ye,logge,dx,dye,dxp,dxm)
          call step(net,xn,xp,ye,logge,dx,dye,dxn,dxp,det,dxdn,dxdp,dyedn,dyedp)
          
          write(99,'(99es13.4e3)') xn,xp,dx,dye,dxn,dxp,det,dxdn,dxdp,dyedn,dyedp
          !write(99,'(99es13.4e3)') xp,xm,dx,dye,dxp,dxm
          
       enddo
    enddo
    close(99)

  end subroutine test_converge


  subroutine step(net,xn,xp,ye,logge,dx,dye,dxn,dxp,det_out,dxdn_out,dxdp_out,dyedn_out,dyedp_out,rcond_out)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    type(nse_network_t),intent(in) :: net
    real(8),intent(in) :: xn,xp,ye
    real(8),intent(in) :: logge(net%n_spec)
    real(8),intent(out) :: dx,dye,dxn,dxp
    real(8),intent(out),optional :: det_out,dxdn_out,dxdp_out,dyedn_out,dyedp_out, rcond_out
    
    real(8) :: logx(net%n_spec)

    real(8) :: det

    real(8) :: u(net%n_spec), logx_max, uq(net%n_spec), logq_max, usum, qsum
    real(8) :: f1, f2, df1dn, df1dp, df2dn, df2dp

    real(8),parameter :: rcond_min = 1d-12
    real(8) :: jnorm1, adjnorm1, rcond

    logx(:) = logge(:) + net%z(:)*xp + net%n(:)*xn

    logx_max = maxval(logx(:))

    u(:) = exp(logx(:) - logx_max)

    usum  = sum(u(:))

    ! residuals
    f1 = logx_max + log(usum)
    
    ! Jacobian
    df1dn = sum(net%n(:)*u(:)) / usum
    df1dp = sum(net%z(:)*u(:)) / usum
    
    logq_max = maxval(logx(:), mask=net%zai(:) > 0d0)
    
    uq(:) = 0d0
    where (net%zai(:) > 0d0)
       uq(:) = net%zai(:)*exp(logx(:) - logq_max)
    end where
    qsum = sum(uq(:))

    f2 = logq_max + log(qsum) - log(ye)

    df2dn = sum(net%n(:)*uq(:))/qsum
    df2dp = sum(net%z(:)*uq(:))/qsum
    
    det = df1dn*df2dp - df1dp*df2dn

    jnorm1 = max(abs(df1dn) + abs(df2dn), abs(df1dp) + abs(df2dp))
    
    adjnorm1 = max(abs(df2dp) + abs(df2dn), abs(df1dp) + abs(df1dn))
    
    if (jnorm1 > 0d0 .and. adjnorm1 > 0d0) then
       rcond = abs(det)/(jnorm1*adjnorm1)
    else
       rcond = 0d0
    endif

    ! if (.not. ieee_is_finite(det) .or. &
    !      .not. ieee_is_finite(rcond)) then
    !    jac_bad = .true.
    !    return
    ! endif
    
    !dxn =-( dx*dyedp-dye*dxdp)/det
    !dxp =-(-dx*dyedn+dye*dxdn)/det
    dx = f1
    dye= f2
    if(present(det_out)) det_out = det
    if(present(dxdn_out)) dxdn_out = df1dn
    if(present(dxdp_out)) dxdp_out = df1dp
    if(present(dyedn_out)) dyedn_out = df2dn
    if(present(dyedp_out)) dyedp_out = df2dp

    if (rcond < rcond_min) then
       ! jac_bad = .true.
       dxn = 0d0
       dxp = 0d0
       if(present(rcond_out)) rcond_out = rcond
       return
    endif

    dxn = (-f1*df2dp + df1dp*f2)/det
    dxp = ( df2dn*f1 - df1dn*f2)/det

  end subroutine step



  ! subroutine step2(xp,xm,ye,logge,dx,dye,dxp,dxm)
  !   real(8),intent(in) :: xm,xp,ye
  !   real(8),intent(in) :: logge(n_spec)
  !   real(8),intent(out) :: dx,dye,dxp,dxm
    
  !   real(8) :: logx(n_spec), x(n_spec)

  !   real(8) :: xsum,yesum,dxdp,dxdm,dyedp,dyedm,det
    
  !   logx(1:n_spec) = logge(1:n_spec) + a(1:n_spec)*xp + (n(1:n_spec)-z(1:n_spec))*xm
    
  !   x(1:n_spec) = 10d0**logx(1:n_spec)
    
  !   xsum = sum(x(1:n_spec))
  !   yesum= sum(zai(1:n_spec)*x(1:n_spec))

  !   dxdp = sum(a(1:n_spec)*x(1:n_spec))
  !   dxdm = sum((n(1:n_spec)-z(1:n_spec))*x(1:n_spec))
    
  !   dyedp = sum(z(1:n_spec)*x(1:n_spec))
  !   dyedm = sum((n(1:n_spec)-z(1:n_spec))*zai(1:n_spec)*x(1:n_spec))
    
  !   dxdp = dxdp/xsum
  !   dxdm = dxdm/xsum

  !   dyedp = dyedp/yesum! - dxdp    
  !   dyedm = dyedm/yesum! - dxdn

  !   dx  = log10(xsum)
  !   dye = log10(yesum/ye)

  !   det = dxdm*dyedp - dxdp*dyedm
  !   if( det == 0d0)then
  !      !write(6,*) "det = 0"
  !      dxp = 0d0
  !      dxm = 0d0
  !   else

  !      dxm =-( dx*dyedp-dye*dxdp)/det
  !      dxp =-(-dx*dyedm+dye*dxdm)/det
  !   endif

  ! end subroutine step2

  subroutine nse_alpha(rho,tem,ye,y_alpha,y_n,y_p)
    use const
    implicit none
    real(8),parameter :: dma = 2.d0*mpmev + 2.d0*mnmev - mamev
    real(8),intent(in) :: rho,ye,tem ! temperature in MeV !
    real(8),intent(out) :: y_alpha,y_n,y_p
    real(8),parameter :: tol = 1.d-15

    integer :: it,itlim
    real(8) :: coef,f,dfdy,dy
    itlim=200


    !coef = 1.6843404717d-35 *rho**3 *tem**(-4.5d0) *exp( min(2.d2,dma/tem) )
    coef = ((2d0*pi*hbc*hbc/mumev)**1.5d0/mu)**3 *rho**3 *tem**(-4.5d0) *exp( min(2.d2,dma/tem) )
    !write(6,*) ((2d0*pi*hbc*hbc/mumev)**1.5d0/mu)**3
    
    if(coef>1.d300) then
       y_alpha = min(ye,1.d0-ye) *0.5d0
       y_n = max(1d-99,1.d0-ye-2.d0*y_alpha)
       y_p = max(1d-99,ye-2.d0*y_alpha)

       !write(6,*) "alpha and one of the nucleons!"
       return
    endif

    y_alpha = min(ye,1.d0-ye) *0.25d0

    do it = 1, itlim
       f    = y_alpha  - 0.5d0*coef *(1.d0-ye-2.d0*y_alpha)**2 *(ye-2.d0*y_alpha)**2
       dfdy = 1.d0 + 2.d0 *coef *(1.d0-ye-2.d0*y_alpha)    *(ye-2.d0*y_alpha)    *(1.d0-4.d0*y_alpha)
       dy   = - f/dfdy
       
       !write(*,*) it,y_alpha,dy,f
       if(isnan(dy))then
          write(*,*) it,y_alpha,dy,f,dfdy,coef
          stop
       endif

       if(abs(dy)<tol) exit
       y_alpha = y_alpha + dy
    enddo

    y_alpha = max(1.d-99,min(ye*0.5d0,y_alpha))

    y_n = max(1d-99,1.d0-ye-2.d0*y_alpha)
    y_p = max(1d-99,ye-2.d0*y_alpha)

    ! if(y_n/y_p<1d0)then
    !    y_n = sqrt(y_alpha/coef)/y_p
    ! elseif(y_p/y_n<1d0)then
    !    y_p = sqrt(y_alpha/coef)/y_n
    ! endif

    ! y_n = max(1d-99,1.d0-ye-2.d0*y_alpha)
    ! y_p = max(1d-99,ye-2.d0*y_alpha)

    return
  end subroutine nse_alpha
  
  subroutine statistic(net, x, mexc_ave, z_heavy, a_heavy, y_heavy, ytot, xsum, yesum)
    type(nse_network_t), intent(in) :: net
    real(8),intent(in) :: x(net%n_spec)
    real(8),intent(out) :: mexc_ave, z_heavy, a_heavy, y_heavy, ytot, xsum, yesum

    integer :: k

    mexc_ave = 0.d0
    do k=1,net%n_spec
       mexc_ave = mexc_ave + net%mexc(k)*x(k)/net%a(k)
    enddo

    ytot = 0.d0
    do k=1,net%n_spec
       ytot = ytot + x(k)/net%a(k)
    enddo

    xsum = 0.d0
    do k=1,net%n_spec
       xsum = xsum + x(k)
    enddo

    yesum = 0.d0
    do k=1,net%n_spec
       yesum = yesum + x(k)/net%a(k)*net%z(k)
    enddo

    z_heavy = 0.d0
    a_heavy = 0.d0
    y_heavy = 0.d0
    do k=1,net%n_spec
       if(net%a(k)>4d0)then
          z_heavy = z_heavy + net%z(k)*x(k)/net%a(k)
          a_heavy = a_heavy + net%a(k)*x(k)/net%a(k)
          y_heavy = y_heavy +          x(k)/net%a(k)
       endif
    enddo
    z_heavy = z_heavy / y_heavy
    a_heavy = a_heavy / y_heavy

  end subroutine statistic


  subroutine statistic_compose(net, rho, x, stat)
    use const, only:emev
    type(nse_network_t),intent(in) :: net
    real(8),intent(in) :: rho, x(net%n_spec)
    type(stat_t),intent(out) :: stat

    integer :: k

    real(8) :: ytot, mexc_ave, yesum, ecoul_ave
    real(8) :: z_heavy, a_heavy, y_heavy

    yesum = 0.d0
    do k=1,net%n_spec
       yesum = yesum + x(k)/net%a(k)*net%z(k)
    enddo

    ! mass-excess per baryon
    ! add m_e*c^2 * Ye to account for the rest-mass of balence electrons.
    mexc_ave = 0.d0
    do k=1,net%n_spec
       mexc_ave = mexc_ave + net%mexc(k)*x(k)/net%a(k)
    enddo
    mexc_ave = mexc_ave + emev*yesum
    stat%mexc = mexc_ave

    ytot = 0.d0
    do k=1,net%n_spec
       ytot = ytot + x(k)/net%a(k)
    enddo
    stat%abar = 1d0/ytot
    
    stat%yn = 0d0
    stat%yp = 0d0
    stat%yh2 = 0d0
    stat%yh3 = 0d0
    stat%yhe3 = 0d0
    stat%yhe4 = 0d0
   
    z_heavy = 0.d0
    a_heavy = 0.d0
    y_heavy = 0.d0
    do k=1,net%n_spec
       ! write(6,*) k, name_nucl(k), x(k), nint(net%a(k)), nint(net%z(k))
       if    (nint(net%a(k))==1.and.nint(net%z(k))==0)then ! n
          stat%yn = x(k)/net%a(k)
       elseif(nint(net%a(k))==1.and.nint(net%z(k))==1)then ! p
          stat%yp = x(k)/net%a(k)
       elseif(nint(net%a(k))==2.and.nint(net%z(k))==1)then ! h2
          stat%yh2 = x(k)/net%a(k)
       elseif(nint(net%a(k))==3.and.nint(net%z(k))==1)then ! h3
          stat%yh3 = x(k)/net%a(k)
       elseif(nint(net%a(k))==3.and.nint(net%z(k))==2)then ! he3
          stat%yhe3 = x(k)/net%a(k)
       elseif(nint(net%a(k))==4.and.nint(net%z(k))==2)then ! he4
          stat%yhe4 = x(k)/net%a(k)
       else
          z_heavy = z_heavy + net%z(k)*x(k)/net%a(k)
          a_heavy = a_heavy + net%a(k)*x(k)/net%a(k)
          y_heavy = y_heavy +          x(k)/net%a(k)
       endif
    enddo
    z_heavy = z_heavy / y_heavy
    a_heavy = a_heavy / y_heavy

    stat%a_n = a_heavy
    stat%z_n = z_heavy
    stat%y_n = y_heavy

    call calc_coulomb_average(net, rho, yesum, x, ecoul_ave)
    
    stat%ecoul = ecoul_ave
  end subroutine statistic_compose


  subroutine output_composition(net, x,temp,rho,ye)
    use module_nuclear_data_winvne
    type(nse_network_t), intent(in) :: net
    real(8),intent(in) :: x(net%n_spec)
    real(8),intent(in) :: temp,rho,ye
    
    real(8) :: mexc_ave, z_heavy, a_heavy, y_heavy, ytot, xsum, yesum
    real(8),allocatable :: xa(:), xz(:), ya(:), yz(:)
    
    integer :: a_max, z_max, ia,iz,k

    a_max = nint(maxval(net%a(:)))
    z_max = nint(maxval(net%z(:)))

    allocate(xa(a_max),ya(a_max),xz(0:z_max),yz(0:z_max))
    
    xa(:)=0d0
    xz(:)=0d0
    ya(:)=0d0
    yz(:)=0d0

    do k=1,net%n_spec
       ia = nint(net%a(k))
       iz = nint(net%z(k))
       xa(ia) = xa(ia) + x(k)
       xz(iz) = xz(iz) + x(k)
       
       ya(ia) = ya(ia) + x(k)/net%a(k)
       yz(iz) = yz(iz) + x(k)/net%a(k)
    enddo

    call statistic(net, x, mexc_ave, z_heavy, a_heavy, y_heavy, ytot, xsum, yesum)

    open(11,file="aabun",status="replace",action="write")
    write(11,'("# T,rho,Ye = ",99es12.4)') temp,rho,ye
    write(11,'("#",99es12.4)') mexc_ave, z_heavy, a_heavy, ytot, xsum
    do ia=1,a_max
       write(11,'(i5,99es15.6e3)') ia, xa(ia), ya(ia)
    enddo
    close(11)

    open(11,file="zabun",status="replace",action="write")
    write(11,'("# T,rho,Ye = ",99es12.4)') temp,rho,ye
    write(11,'("#",99es12.4)') mexc_ave, z_heavy, a_heavy, ytot, xsum
    do iz=0,z_max
       write(11,'(i5,99es15.6e3)') iz, xz(iz), yz(iz)
    enddo
    close(11)


    open(11,file="abun",status="replace",action="write")
    write(11,'("# T,rho,Ye = ",99es12.4)') temp,rho,ye
    write(11,'("#",99es12.4)') mexc_ave, z_heavy, a_heavy, ytot, xsum
    do k=1,net%n_spec
       write(11,'(3i5,99es15.6e3)') nint(net%n(k)),nint(net%z(k)),nint(net%a(k)),x(k),x(k)/net%a(k)
    enddo
    ! do ia=1,a_max
    !    write(11,*)
    !    !do iz=0,min(z_max,ia-1)
    !    do iz=0,z_max
    !       if(jnuc_winvne(ia,iz)==0)then
    !          x_dummy = 0d0
    !          y_dummy = 0d0
    !       else
    !          x_dummy = x(jnuc_winvne(ia,iz))
    !          y_dummy = x(jnuc_winvne(ia,iz))/dble(ia)
    !       endif
    !       write(11,'(2i5,99es15.6e3)') iz,ia-iz, x_dummy,y_dummy
    !    enddo
    ! enddo
    close(11)


  end subroutine output_composition
  
  subroutine output_nse_full(net,rho,temp,ye,xnse, output_filename)
    type(nse_network_t),intent(in) :: net
    real(8),intent(in) :: rho,temp,ye
    real(8),intent(in) :: xnse(net%n_spec)
    character(*),intent(in) :: output_filename
    integer :: unit
    integer :: i
    real(8) :: g(net%n_spec)
    real(8) :: t9
    
    t9 = temp/1d9
    call calc_ptf_nse(net,t9,g)
    open(newunit=unit, file=output_filename, status="replace", action="write")
    write(unit, '("#",a,es15.7,a,es15.7,a,es15.7)') "rho=",rho,"T=",temp,"ye=",ye
    write(unit, '("#",99a20)') "index", "name", "A", "Z", "N", "Xi", "Yi", "gi"
    do i=1,net%n_spec
       write(unit, '(" ",i20,a20,3i20,3es20.10e3)') i, net%name_nucl(i), nint(net%a(i)), nint(net%z(i)), nint(net%n(i)), xnse(i), xnse(i)/net%a(i), g(i)
    end do
    close(unit)
  end subroutine output_nse_full
  
  subroutine calc_coulomb_average(net, rho, ye, x, ecoul_ave)
    type(nse_network_t),intent(in) :: net
    real(8),intent(in)  :: rho, ye
    real(8),intent(in)  :: x(net%n_spec)
    real(8),intent(out) :: ecoul_ave
    
    real(8) :: fcoul(net%n_spec)
    
    call calc_coulomb_HS(net, rho, ye, fcoul)
    
    ecoul_ave = sum(x(:)/net%a(:) * fcoul(:))
    
  end subroutine calc_coulomb_average

  subroutine resolve_network_stat_weights(net, iwinvne, irauscher)
    
    type(nse_network_t), intent(inout) :: net
    integer, intent(in) :: iwinvne(net%n_spec)
    integer, intent(in) :: irauscher(net%n_spec)

    integer :: i

    allocate(net%stat_weight(net%n_spec))

    do i = 1, net%n_spec

       call resolve_stat_weight_ref( &
            net%stat_weight_policy, &
            iwinvne(i), &
            irauscher(i), &
            net%stat_weight(i))

       if (net%stat_weight(i)%source == STAT_WEIGHT_NONE) then
          write(*,*) "ERROR: no statistical-weight source for ", &
               net%name_nucl(i)
          error stop
       endif

    enddo

  end subroutine resolve_network_stat_weights

end module module_nse
