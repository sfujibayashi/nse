module module_nse
  implicit none

  private
  public :: nse_init_four,calc_nse,test_converge,nse_init_reaclib,output_composition,statistic, two_nuclei_approx

  integer :: n_spec
  real(8),allocatable :: mexc(:), a(:), z(:), n(:), g(:), zai(:)
  character(5),allocatable :: name_nucl(:)
  
  logical :: use_reaclib
  
contains
  
  subroutine nse_init_four(n_spec_out)
    use const,only:mnmev,mpmev,mamev,mumev,memev
    integer, intent(out) :: n_spec_out
    real(8),parameter :: mexc_56ni_mev = -53.907539d0

    use_reaclib = .false.

    n_spec = 4
    allocate(mexc(n_spec), a(n_spec), z(n_spec), n(n_spec), g(n_spec),zai(n_spec))
    ! n
    a(1) = 1d0; z(1) = 0d0; n(1) = 1d0; g(1) = 2d0; mexc(1) = mnmev-a(1)*mumev
    ! p
    a(2) = 1d0; z(2) = 1d0; n(2) = 0d0; g(2) = 2d0; mexc(2) = mpmev-a(2)*mumev
    ! alpha
    a(3) = 4d0; z(3) = 2d0; n(3) = 2d0; g(3) = 1d0; mexc(3) = mamev-a(3)*mumev
    ! 56Ni
    a(4) =56d0; z(4) =28d0; n(4) =28d0; g(4) = 1d0; mexc(4) = mexc_56ni_mev-z(4)*memev

    zai(1:n_spec) = z(1:n_spec)/a(1:n_spec)
    
    n_spec_out = n_spec
    
  end subroutine nse_init_four

  subroutine nse_init_reaclib(n_spec_out)
    use module_ptf_reaclib
    use const,only:memev
    integer, intent(out) :: n_spec_out
    integer :: k

    use_reaclib = .true.

    n_spec = nct_reaclib
    allocate(name_nucl(n_spec))
    allocate(mexc(n_spec), a(n_spec), z(n_spec), n(n_spec), g(n_spec),zai(n_spec))

    do k=1,nct_reaclib
       a(k) = ams_reaclib(k)
       z(k) = dble(npt_reaclib(k))
       n(k) = dble(nnt_reaclib(k))
       mexc(k) = exc_reaclib(k) - z(k)*memev
       name_nucl(k) = name_reaclib(k)
    enddo

    ! block 
    !   use const,only : mumev
    !   real(8) :: mexcpb = 0d0
    !   integer :: k_min
    !   k_min = 0
    !   do k=1,nct_reaclib
    !      if(mexcpb> mexc(k)/a(k))then
    !         mexcpb = mexc(k)/a(k)
    !         k_min = k
    !      endif
         
    !   enddo
    !   write(6,*) mexcpb, mexcpb/mumev, a(k_min),z(k_min)
    !   stop
    ! end block

    zai(1:n_spec) = z(1:n_spec)/a(1:n_spec)
    
    n_spec_out = n_spec
    
  end subroutine nse_init_reaclib

  subroutine calc_coulomb(rho,ye,fcoul)
    use const, only : qe, mu, pi, mev2erg
    real(8),intent(in) :: rho,ye
    real(8),intent(out) :: fcoul(n_spec)

    integer :: k
    real(8) :: ne, v_n, v_c, u
    real(8),parameter :: n0 = 0.16d0*1d39
    
    fcoul(:) = 0d0
    
    ne = ye*rho/mu
    do k=1,n_spec
       if( z(k)>0.d0)then
          v_n = a(k)/n0
          v_c = z(k)/ne
          u   = v_n/v_c
          
          fcoul(k) = (3d0/5d0)*(4d0*pi/3d0)**(-1d0/3d0) * qe**2 * n0**2 * (z(k)/a(k))**2 * (v_n)**(5d0/3d0) &
               * (-3d0/2d0*u**(1d0/3d0) + 1d0/2d0*u) / mev2erg
       endif
    enddo
    
  end subroutine calc_coulomb

  function Ecoul_HS10_eq6(z, a, n0, ne) result(ecoul)
    use const, only : pi, fine
    real(8),intent(in) :: z,a,n0,ne
    real(8) :: ecoul
    real(8) :: x,r

    x = (ne/n0 * a/z)**(1d0/3d0)
    r = ((3d0*a)/(4d0*pi*n0))**(1d0/3d0)
    ecoul = -3d0/5d0 * z*z*fine/r*(3d0/2d0*x - 1d0/2d0*x*x*x) 
    
  end function Ecoul_HS10_eq6

  subroutine calc_ptf_HS10(t9,g)
    real(8),intent(in) :: t9
    real(8),intent(out) :: g(n_spec)
    
    integer :: k

    do k=1,n_spec
       
    end do
    
  end subroutine calc_ptf_HS10

  ! function excited_HS10(temp, a) result g
  !   real(8),intent(in) :: temp, a
  !   real(8),parameter :: c1=0.2, c2=0.8
    
  ! end function excited_HS10
  
  subroutine calc_nse(rho,temp,ye,itrlim,tol,xnse,nsefail,use_TNAguess,xn_history,xp_history,itr_out,err_out,xn_guess,xp_guess,xn_out,xp_out)
    use const,only : mu,kerg,pi,hbar,mev2erg
    use module_ptf_reaclib
    real(8),intent(in) :: rho,temp,ye
    integer,intent(in) :: itrlim
    real(8),intent(in) :: tol
    real(8),intent(out) :: xnse(n_spec)
    logical,intent(out) :: nsefail
    logical,intent(in) :: use_TNAguess
    real(8),intent(out),optional :: xn_history(0:itrlim),xp_history(0:itrlim)
    integer,intent(out),optional :: itr_out
    real(8),intent(out),optional :: err_out
    real(8),intent(in),optional :: xn_guess,xp_guess
    real(8),intent(out),optional :: xn_out,xp_out
    
    real(8) :: logrho0
    real(8) :: logge(n_spec), logx(n_spec), x(n_spec), fcoul(n_spec)
    
    real(8) :: xp,xn

    integer :: itr
    !integer,parameter :: itrlim=50
    !    real(8),parameter :: tol = 1d-13
    
    real(8) :: dxdp,dxdn,dyedp,dyedn,dx,dye,det,dxn,dxp,dl,fac

    real(8) :: t9

    real(8),parameter :: n0 = 0.16d0*1d39
    integer :: i_spec
    real(8) :: nb
    
    ! log(rho0/rho)
    logrho0 = 2.5d0*log(mu) + 1.5d0*log(kerg*temp) - 1.5d0*log(2d0*pi) - 3d0*log(hbar) - log(rho)
    
    ! partition function may be calculated here
    t9 = temp/1d9
    if(use_reaclib)call calc_ptf(t9,g)

    if(use_reaclib)then
       call calc_coulomb(rho,ye,fcoul)
    else
       fcoul(:) = 0d0
    endif
    !
    logge(1:n_spec) = log(g(1:n_spec)) + 2.5d0*log(a(1:n_spec)) + logrho0 - mexc(1:n_spec)*mev2erg/(kerg*temp) &
         - fcoul(1:n_spec)*mev2erg/(kerg*temp)
    
    ! nb = rho/mu
    ! do i_spec = 1,n_spec
    !    if(a(i_spec) > 1)then
    !       logge(i_spec) = logge(i_spec) + log10(max(1d-99,1d0-nb/n0))
    !    endif
    ! enddo
    
    ! write(6,'(99es12.4)') logrho0, log10(mu*(mu*kerg*temp/(2d0*pi*hbar*hbar))**1.5d0/rho)
    ! write(6,'(99es12.4)') logge(1:n_spec)
    ! write(6,'(99es12.4)') mexc(1:n_spec)
    ! write(6,'(99es12.4)') - mexc(1:n_spec)*mev2erg/(kerg*temp)/log(10d0)
    ! stop
    
    if(present(xn_history)) xn_history(:) = 0d0
    if(present(xp_history)) xp_history(:) = 0d0


    ! use two-nuclei approx for initial guess
    if(use_TNAguess)then
       block
         integer :: k1,k2
         real(8) :: eta01ex, eta02ex
         real(8) :: z1,z2,a1,a2,x1,x2,g1,g2,mex1,mex2,n1,n2

         if(ye/=0.5d0)then
            call two_nuclei_approx_index(ye, k1, k2)
         else
            k1 = jnuc_reaclib(1,1)
            k2 = jnuc_reaclib(56,26)
         endif

         z1 = z(k1)
         z2 = z(k2)
         a1 = a(k1)
         a2 = a(k2)
         n1 = a1-z1
         n2 = a2-z2
         x2 = (z1/a1 - ye)/(z1/a1 - z2/a2)
         x1 = (1d0 - x2)
         
         g1 = g(k1)
         g2 = g(k2)
         mex1 = mexc(k1)*mev2erg
         mex2 = mexc(k2)*mev2erg
         
         ! (mu_1 - m_1 c^2 + mexc_1*c^2)/kT
         eta01ex = -logrho0 + log(x1) - log(g1) - 2.5d0*log(a1) + (mexc(k1) + fcoul(k1))*mev2erg/(kerg*temp)
         ! (mu_2 - m_2 c^2 + mexc_2*c^2)/kT
         eta02ex = -logrho0 + log(x2) - log(g2) - 2.5d0*log(a2) + (mexc(k2) + fcoul(k2))*mev2erg/(kerg*temp)
         
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
             call step(xn,xp,ye,logge,dx,dye,dxn,dxp,det,dxdn,dxdp,dyedn,dyedp)
             logx(1:n_spec) = logge(1:n_spec) + z(1:n_spec)*xp + n(1:n_spec)*xn
             !write(6,'(99es12.4)') xn,xp,det,maxval(logx(:)),minval(logx(:))
             if( abs(det)>0d0 .and. maxval(logx(:))<3d2*log(10d0) .and. dx<0d0 .and.dye<0d0)then
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

       call step(xn,xp,ye,logge,dx,dye,dxn,dxp,det,dxdn,dxdp,dyedn,dyedp)
       
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

    logx(1:n_spec) = logge(1:n_spec) + z(1:n_spec)*xp + n(1:n_spec)*xn 
    logx(1:n_spec) = max(-3d2*log(10d0),min(3d2*log(10d0),logx(1:n_spec)))
    
    x(1:n_spec) = exp(logx(1:n_spec))
    !write(6,*) itr,xn,xp
    xnse(:) = x(:)

    if(present(xn_out)) xn_out = xn
    if(present(xp_out)) xp_out = xp

  end subroutine calc_nse

  subroutine two_nuclei_approx(ye,xnse)
    real(8),intent(in) :: ye
    real(8),intent(out) :: xnse(n_spec)
    
    real(8) :: y1,y2, z1,z2,a1,a2
    integer :: k1,k2

    call two_nuclei_approx_index(ye, k1, k2)
    
    ! real(8) :: y1,y2, z1,z2,a1,a2, mexc1, mexc2, f, f_min


    ! integer :: k1_min, k2_min

    ! k1_min = 0
    ! k2_min = 0
    ! f_min = 1d99
    ! do k1=1,n_spec
    !    do k2=1,k1-1
    !       z1 = z(k1)
    !       z2 = z(k2)
    !       a1 = a(k1)
    !       a2 = a(k2)
    !       if( (z1/a1 - ye)*(z1/a1 - z2/a2) > 0d0 )then
    !          y2 = (z1/a1 - ye)/(z1/a1 - z2/a2)/a2
    !          y1 = (1d0 - a2*y2)/a1
             
    !          mexc1 = mexc(k1)
    !          mexc2 = mexc(k2)
             
    !          f = mexc1*y1 + mexc2*y2

    !          !write(6,*) k1,k2,y1,y2,f
    !          if(y1>0d0 .and. f < f_min)then
    !             f_min = f
    !             k1_min = k1
    !             k2_min = k2
                
    !             !write(6,'(2i5,99es12.4)') k1_min,k2_min,f_min,mexc1,mexc2, y1,y2
    !             !stop
    !          endif
    !       endif
    !    enddo
    ! enddo
    
    ! k1 = k1_min
    ! k2 = k2_min

    z1 = z(k1)
    z2 = z(k2)
    a1 = a(k1)
    a2 = a(k2)
    y2 = (z1/a1 - ye)/(z1/a1 - z2/a2)/a2
    y1 = (1d0 - a2*y2)/a1
    
    xnse(:) = 1d-300
    xnse(k1) = a1*y1
    xnse(k2) = a2*y2
    
  end subroutine two_nuclei_approx


  subroutine two_nuclei_approx_index(ye,k1_min,k2_min)
    real(8),intent(in) :: ye
    integer,intent(out) :: k1_min,k2_min
    real(8) :: y1,y2, z1,z2,a1,a2, mexc1, mexc2, f, f_min

    integer :: k1, k2

    k1_min = 0
    k2_min = 0
    f_min = 1d99
    do k1=1,n_spec
       do k2=1,k1-1
          z1 = z(k1)
          z2 = z(k2)
          a1 = a(k1)
          a2 = a(k2)
          if( (z1/a1 - ye)*(z1/a1 - z2/a2) > 0d0 )then
             y2 = (z1/a1 - ye)/(z1/a1 - z2/a2)/a2
             y1 = (1d0 - a2*y2)/a1
             
             mexc1 = mexc(k1)
             mexc2 = mexc(k2)
             
             f = mexc1*y1 + mexc2*y2

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

  subroutine test_converge(rho,temp,ye,use_TNAguess)
    use const,only : mu,kerg,pi,hbar,mev2erg
    use module_ptf_reaclib

    real(8),intent(in) :: rho,temp,ye
    logical,intent(in) :: use_TNAguess
    real(8) :: logrho0
    real(8) :: logge(n_spec), fcoul(n_spec)
    
    real(8) :: xp,xn

    real(8) :: dxdp,dxdn,dyedp,dyedn,dx,dye,det,dxn,dxp

    integer :: in,ip,nn,np
    real(8) :: xn_min,xn_max,xp_min,xp_max

    integer :: itr,itr_out
    real(8) :: tol = 1d-12
    integer,parameter :: itrlim=1000
    real(8) :: xn_history(0:itrlim),xp_history(0:itrlim),xnse(n_spec)

    ! real(8) :: xm,dxm,xm_min,xm_max
    logical :: nsefail
    real(8) :: t9

    call calc_nse(rho,temp,ye,itrlim,tol,xnse,nsefail,use_TNAguess,xn_history,xp_history,itr_out)
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
    if(use_reaclib)call calc_ptf(t9,g)

    if(use_reaclib)then
       call calc_coulomb(rho,ye,fcoul)
    else
       fcoul(:) = 0d0
    endif
    !
    logge(1:n_spec) = log(g(1:n_spec)) + 2.5d0*log(a(1:n_spec)) + logrho0 - mexc(1:n_spec)*mev2erg/(kerg*temp) &
         - fcoul(1:n_spec)*mev2erg/(kerg*temp)

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
          call step(xn,xp,ye,logge,dx,dye,dxn,dxp,det,dxdn,dxdp,dyedn,dyedp)
          
          write(99,'(99es13.4e3)') xn,xp,dx,dye,dxn,dxp,det,dxdn,dxdp,dyedn,dyedp
          !write(99,'(99es13.4e3)') xp,xm,dx,dye,dxp,dxm
          
       enddo
    enddo
    close(99)

  end subroutine test_converge


  subroutine step(xn,xp,ye,logge,dx,dye,dxn,dxp,det_out,dxdn_out,dxdp_out,dyedn_out,dyedp_out)
    real(8),intent(in) :: xn,xp,ye
    real(8),intent(in) :: logge(n_spec)
    real(8),intent(out) :: dx,dye,dxn,dxp
    real(8),intent(out),optional :: det_out,dxdn_out,dxdp_out,dyedn_out,dyedp_out
    
    real(8) :: logx(n_spec), x(n_spec)

    real(8) :: xsum,yesum
    real(8) :: det,dxdn,dxdp,dyedn,dyedp

    real(8) :: u(n_spec), logx_max, logx_sum, usum, w(n_spec), qbar
    real(8) :: f1, f2, nbar, zbar, qnbar, qzbar, df1dn, df1dp, df2dn, df2dp
    ! real(8) :: ave_a,a11,a12,a21,a22

    logx(1:n_spec) = logge(1:n_spec) + z(1:n_spec)*xp + n(1:n_spec)*xn

    logx_max = maxval(logx(:))
    
    u(:) = exp(logx(:) - logx_max)
    usum = sum(u(:))

    w(:) = u(:)/usum

    logx_sum = logx_max + log(usum)

    x(1:n_spec) = exp(logx(1:n_spec))
    
    qbar = sum(zai(:)*w(:))

    ! residuals
    f1 = logx_sum
    f2 = qbar - ye

    ! weighted moments
    nbar  = sum(n(:)*w(:))
    zbar  = sum(z(:)*w(:))
    qnbar = sum(zai(:)*n(:)*w(:))
    qzbar = sum(zai(:)*z(:)*w(:))

    ! Jacobian
    df1dn = nbar
    df1dp = zbar

    df2dn = qnbar - qbar*nbar
    df2dp = qzbar - qbar*zbar
    
    det = df1dn*df2dp - df1dp*df2dn
    if( det == 0d0)then
       !write(6,*) "det = 0",dxdn*dyedp, dxdp*dyedn, dxdn,dxdp,dyedn,dyedp
       !det = 1d0
       dxn = (-f1*df2dp + df1dp*f2)
       dxp = ( df2dn*f1 - df1dn*f2)
       ! ave_a = 2d0
       ! a11 = - (x(1)+2d0*x(2))/x(3)/ave_a
       ! a12 = - (x(2)+2d0*x(1))/x(3)/ave_a
       ! a21 = - (     4d0*x(2))/x(3)/ave_a
       ! a22 = - (     2d0*x(2))/x(3)/ave_a
       ! det = ave_a**2*(a11+a22+a11*a22 - (a12+a21+a12*a21))
       
       !write(6,*) det
       !stop
    else
       dxn = (-f1*df2dp + df1dp*f2)/det
       dxp = ( df2dn*f1 - df1dn*f2)/det
    endif

    !dxn =-( dx*dyedp-dye*dxdp)/det
    !dxp =-(-dx*dyedn+dye*dxdn)/det
    dx = f1
    dye= f2
    if(present(det_out)) det_out = det
    if(present(dxdn_out)) dxdn_out = df1dn
    if(present(dxdp_out)) dxdp_out = df1dp
    if(present(dyedn_out)) dyedn_out = df2dn
    if(present(dyedp_out)) dyedp_out = df2dp
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
  
  subroutine statistic(x,mexc_ave, z_heavy, a_heavy, y_heavy, ytot, xsum, yesum)
    real(8),intent(in) :: x(n_spec)
    real(8),intent(out) :: mexc_ave, z_heavy, a_heavy, y_heavy, ytot, xsum, yesum

    integer :: k

    mexc_ave = 0.d0
    do k=1,n_spec
       mexc_ave = mexc_ave + mexc(k)*x(k)/a(k)
    enddo

    ytot = 0.d0
    do k=1,n_spec
       ytot = ytot + x(k)/a(k)
    enddo

    xsum = 0.d0
    do k=1,n_spec
       xsum = xsum + x(k)
    enddo

    yesum = 0.d0
    do k=1,n_spec
       yesum = yesum + x(k)/a(k)*z(k)
    enddo

    z_heavy = 0.d0
    a_heavy = 0.d0
    y_heavy = 0.d0
    do k=1,n_spec
       if(a(k)>4d0)then
          z_heavy = z_heavy + z(k)*x(k)/a(k)
          a_heavy = a_heavy + a(k)*x(k)/a(k)
          y_heavy = y_heavy +      x(k)/a(k)
       endif
    enddo
    z_heavy = z_heavy / y_heavy
    a_heavy = a_heavy / y_heavy

  end subroutine statistic

  subroutine output_composition(x,temp,rho,ye)
    use module_ptf_reaclib
    real(8),intent(in) :: x(n_spec)
    real(8),intent(in) :: temp,rho,ye
    
    real(8) :: mexc_ave, z_heavy, a_heavy, y_heavy, ytot, xsum, yesum
    real(8),allocatable :: xa(:), xz(:), ya(:), yz(:)
    
    integer :: a_max, z_max, ia,iz,k

    a_max = nint(maxval(a(:)))
    z_max = nint(maxval(z(:)))

    allocate(xa(a_max),ya(a_max),xz(0:z_max),yz(0:z_max))
    
    xa(:)=0d0
    xz(:)=0d0
    ya(:)=0d0
    yz(:)=0d0

    do k=1,n_spec
       ia = nint(a(k))
       iz = nint(z(k))
       xa(ia) = xa(ia) + x(k)
       xz(iz) = xz(iz) + x(k)
       
       ya(ia) = ya(ia) + x(k)/a(k)
       yz(iz) = yz(iz) + x(k)/a(k)
    enddo

    call statistic(x,mexc_ave, z_heavy, a_heavy, y_heavy, ytot, xsum, yesum)

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
    do k=1,n_spec
       write(11,'(3i5,99es15.6e3)') nint(n(k)),nint(z(k)),nint(a(k)),x(k),x(k)/a(k)
    enddo
    ! do ia=1,a_max
    !    write(11,*)
    !    !do iz=0,min(z_max,ia-1)
    !    do iz=0,z_max
    !       if(jnuc_reaclib(ia,iz)==0)then
    !          x_dummy = 0d0
    !          y_dummy = 0d0
    !       else
    !          x_dummy = x(jnuc_reaclib(ia,iz))
    !          y_dummy = x(jnuc_reaclib(ia,iz))/dble(ia)
    !       endif
    !       write(11,'(2i5,99es15.6e3)') iz,ia-iz, x_dummy,y_dummy
    !    enddo
    ! enddo
    close(11)


  end subroutine output_composition
  

end module module_nse
