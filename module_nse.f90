module module_nse
  implicit none

  private
  public :: nse_init,calc_nse,test_converge,nse_alpha

  integer :: n_spec
  real(8),allocatable :: mexc(:), a(:), z(:), n(:), g(:), zai(:)
  
contains
  
  subroutine nse_init(n_spec_out)
    use const,only:mnmev,mpmev,mamev,mumev
    integer, intent(out) :: n_spec_out
    real(8),parameter :: m56ni_mev = -53.907539d0

    n_spec = 3
    allocate(mexc(n_spec), a(n_spec), z(n_spec), n(n_spec), g(n_spec),zai(n_spec))
    ! n
    a(1) = 1d0; z(1) = 0d0; n(1) = 1d0; g(1) = 2d0; mexc(1) = mnmev-a(1)*mumev
    ! p
    a(2) = 1d0; z(2) = 1d0; n(2) = 0d0; g(2) = 2d0; mexc(2) = mnmev-a(2)*mumev
    ! alpha
    a(3) = 4d0; z(3) = 2d0; n(3) = 2d0; g(3) = 1d0; mexc(3) = mamev-a(3)*mumev
    ! 56Ni
    !a(4) =56d0; z(4) =28d0; n(4) =28d0; g(4) = 1d0; mexc(4) = m56ni_mev-a(4)*mumev
    
    zai(1:n_spec) = z(1:n_spec)/a(1:n_spec)
    
    
    n_spec_out = n_spec
    
  end subroutine nse_init
  
  subroutine calc_nse(rho,temp,ye,itrlim,tol,xnse,xn_history,xp_history,itr_out)
    use const,only : mu,kerg,pi,hbar,mev2erg
    real(8),intent(in) :: rho,temp,ye
    integer,intent(in) :: itrlim
    real(8),intent(in) :: tol
    real(8),intent(out) :: xnse(n_spec)
    real(8),intent(out),optional :: xn_history(0:itrlim),xp_history(0:itrlim)
    integer,intent(out),optional :: itr_out
    
    real(8) :: logrho0
    real(8) :: logge(n_spec), logx(n_spec), x(n_spec)
    
    real(8) :: xp,xn

    integer :: itr
    !integer,parameter :: itrlim=50
    !    real(8),parameter :: tol = 1d-13

    real(8) :: xsum,yesum,dxdp,dxdn,dyedp,dyedn,dx,dye,det,dxn,dxp,dl,fac

    ! log10(rho0/rho)
    logrho0 = 2.5d0*log10(mu) + 1.5d0*log10(kerg*temp) - 1.5d0*log10(2d0*pi) - 3d0*log10(hbar) - log10(rho)
    
    ! partition function may be calculated here
    
    !
    logge(1:n_spec) = log10(g(1:n_spec)) + 2.5d0*log10(a(1:n_spec)) + logrho0 - mexc(1:n_spec)*mev2erg/(kerg*temp)/log(10d0)  
    
    ! write(6,'(99es12.4)') logrho0, log10(mu*(mu*kerg*temp/(2d0*pi*hbar*hbar))**1.5d0/rho)
    write(6,'(99es12.4)') logge(1:n_spec)
    ! stop
    if(present(xn_history)) xn_history(:) = 0d0
    if(present(xp_history)) xp_history(:) = 0d0

    xn = -2d0 - logge(1)
    xp = -2d0 - logge(2)
    do itr=1,10000
       call step(xn,xp,ye,logge,dx,dye,dxn,dxp,det,dxdn,dxdp,dyedn,dyedp)
       logx(1:n_spec) = logge(1:n_spec) + z(1:n_spec)*xp + n(1:n_spec)*xn
       !write(6,'(99es12.4)') xn,xp,det,logx(:)
       if( abs(det)>0d0 .and. maxval(logx(:))<3d2 .and. dx<0d0 .and.dye<0d0)then
          exit
       else
          xn = xn - 1d0
          xp = xp - 1d0
       endif
       
    enddo

    !write(6,'(99es12.4)') xn,xp,det,logx(:)
    call step(xn,xp,ye,logge,dx,dye,dxn,dxp,det,dxdn,dxdp,dyedn,dyedp)
    !write(6,'(99es12.4)') xn,xp,ye,logge,dx,dye,dxn,dxp,det,dxdn,dxdp,dyedn,dyedp

    if(abs(det)/min(abs(dxdn),abs(dxdp),abs(dyedn),abs(dyedp))<1d-15)then
       xnse(3) = min(ye,1d0-ye)*2d0
       xnse(1) = max(1d-99, 1d0-ye - 0.5d0*xnse(3))
       xnse(2) = max(1d-99, ye     - 0.5d0*xnse(3))
       if( abs(ye-0.5d0)<1d-16 )then
          xnse(1) = 1d-99
          xnse(2) = 1d-99
          xnse(3) = 1d0
       endif
       if(present(itr_out)) itr_out = 0
       if(present(xn_history)) xn_history(0) = 0d0
       if(present(xp_history)) xp_history(0) = 0d0
       return
    endif
    
    if(present(xn_history)) xn_history(0) = xn
    if(present(xp_history)) xp_history(0) = xp
    if(present(itr_out))itr_out = 0
    
    do itr=1,itrlim

       if(present(xn_history))xn_history(itr) = xn
       if(present(xp_history))xp_history(itr) = xp
       if(present(itr_out))itr_out = itr

       call step(xn,xp,ye,logge,dx,dye,dxn,dxp,det,dxdn,dxdp,dyedn,dyedp)
       
       write(6,'(i5,99es15.7)') itr,xn,xp,dx,dye,dxn,dxp,det,dxdn,dxdp,dyedn,dyedp
       
       !if( (abs(dx)<tol .and. abs(dye) < tol) .or. (abs(dxn/xn)<tol .and. abs(dxp/xp)<tol) ) exit
       if( abs(dx) < tol .and. abs(dye) < tol ) exit
       
       if( abs(det)/min(abs(dxdn),abs(dxdp),abs(dyedn),abs(dyedp))<1d-15 .or. logge(1) + xn < -3d2 .or. logge(2) + xp < -3d2 )then
          xnse(3) = min(ye,1d0-ye)*2d0
          xnse(1) = max(1d-99, 1d0-ye - 0.5d0*xnse(3))
          xnse(2) = max(1d-99, ye     - 0.5d0*xnse(3))
          if( abs(ye-0.5d0)<1d-16 )then
             xnse(1) = 1d-99
             xnse(2) = 1d-99
             xnse(3) = 1d0
          endif
          return
       endif
       
       dl = sqrt(dxn*dxn+dxp*dxp)
       fac = 1d0
       if(dl>2d0)fac = 2d0/dl
       !if(max(abs(dxn/xn),abs(dxp/xp)) > 0.5d0) fac = 0.5d0/max(abs(dxn/xn),abs(dxp/xp))
       !write(6,*) xn,xp

       xn = xn + dxn*fac
       xp = xp + dxp*fac

       !itr_out = itr
    enddo


    logx(1:n_spec) = logge(1:n_spec) + z(1:n_spec)*xp + n(1:n_spec)*xn 
    logx(1:n_spec) = max(-3d2,min(3d2,logx(1:n_spec)))

    x(1:n_spec) = 10d0**logx(1:n_spec)
    ! write(6,*) xn,xp
    write(6,*) x(:)
    xnse(:) = x(:)

  end subroutine calc_nse

  subroutine test_converge(rho,temp,ye)
    use const,only : mu,kerg,pi,hbar,mev2erg

    real(8),intent(in) :: rho,temp,ye

    real(8) :: logrho0
    real(8) :: logge(n_spec), logx(n_spec), x(n_spec)
    
    real(8) :: xp,xn

    real(8) :: xsum,yesum,dxdp,dxdn,dyedp,dyedn,dx,dye,det,dxn,dxp

    integer :: in,ip,nn,np
    real(8) :: xn_min,xn_max,xp_min,xp_max

    integer :: itr,itr_out
    real(8) :: tol = 1d-15
    integer,parameter :: itrlim=200
    real(8) :: xn_history(0:itrlim),xp_history(0:itrlim),xnse(n_spec)

    real(8) :: xm,dxm,xm_min,xm_max

    call calc_nse(rho,temp,ye,itrlim,tol,xnse,xn_history,xp_history,itr_out)
    write(6,*) xnse(:)
    write(6,*) itr_out

    do itr = 0,itr_out
       write(98,*) itr,xn_history(itr),xp_history(itr)
    enddo
    
    logrho0 = 2.5d0*log10(mu) + 1.5d0*log10(kerg*temp) - 1.5d0*log10(2d0*pi) - 3d0*log10(hbar) - log10(rho)

    logge(1:n_spec) = log10(g(1:n_spec)) + 2.5d0*log10(a(1:n_spec)) + logrho0 - mexc(1:n_spec)*mev2erg/(kerg*temp)/log(10d0)  

    nn=200
    np=200
    xn_min = -30d0
    xn_max =  30d0
    xp_min = -30d0
    xp_max =  30d0

    if(xn_min>xn_history(itr_out)) xn_min = xn_history(itr_out)-5d0
    if(xp_min>xp_history(itr_out)) xp_min = xp_history(itr_out)-5d0
    if(xn_max<xn_history(itr_out)) xn_max = xn_history(itr_out)+5d0
    if(xp_max<xp_history(itr_out)) xp_max = xp_history(itr_out)+5d0
    
    !xm_min = -20d0
    !xm_max = 20d0
    
    write(99,'("#",99es12.4)') rho,temp,ye
    write(99,'("#",99es12.4)') xnse(1:3)
    
    do ip=1,np
       write(99,*)
       do in=1,nn

          xn = xn_min + (xn_max-xn_min)*dble(in-1)/dble(nn-1)
          xp = xp_min + (xp_max-xp_min)*dble(ip-1)/dble(np-1)

          !write(6,*) in,ip,xn,xp

          !xm = xm_min + (xm_max-xm_min)*dble(in-1)/dble(nn-1)

          !call step2(xp,xm,ye,logge,dx,dye,dxp,dxm)
          call step(xn,xp,ye,logge,dx,dye,dxn,dxp,det,dxdn,dxdp,dyedn,dyedp)
          
          !xn = 12.778382506817715d0
          !xp = -7.2831081174406629 
          
          ! logx(1:n_spec) = logge(1:n_spec) + z(1:n_spec)*xp + n(1:n_spec)*xn 
          
          ! x(1:n_spec) = 10d0**logx(1:n_spec)
          
          ! !write(6,*) x(:)

          ! xsum = sum(x(1:n_spec))
          ! yesum= sum(zai(1:n_spec)*x(1:n_spec))

          ! dxdn = sum(n(1:n_spec)*x(1:n_spec))
          ! dxdp = sum(z(1:n_spec)*x(1:n_spec))
          
          ! dyedn = sum(n(1:n_spec)*zai(1:n_spec)*x(1:n_spec))
          ! dyedp = sum(z(1:n_spec)*zai(1:n_spec)*x(1:n_spec))
          
          ! dxdn = dxdn/xsum
          ! dxdp = dxdp/xsum
          
          ! dyedn = dyedn/yesum! - dxdn
          ! dyedp = dyedp/yesum! - dxdp
          
          ! dx  = log10(xsum)
          ! dye = log10(yesum/ye)

          ! det = dxdn*dyedp - dxdp*dyedn
          ! !write(6,*) det
          ! !dxn =-( dx*dyedp-dye*dxdp)/det
          ! !dxp =-(-dx*dyedn+dye*dxdn)/det
          ! if( det == 0d0)then
          !    dxn = 0d0
          !    dxp = 0d0
          ! else
          !    dxn =-( dx*dyedp-dye*dxdp)/det
          !    dxp =-(-dx*dyedn+dye*dxdn)/det
          ! endif
          
          !write(6,'(99es13.4e3)') xsum,yesum,xn,xp,dx,dye,dxn,dxp
          !stop

          write(99,'(99es13.4e3)') xn,xp,dx,dye,dxn,dxp,det,dxdn,dxdp,dyedn,dyedp
          !write(99,'(99es13.4e3)') xp,xm,dx,dye,dxp,dxm
          
       enddo
    enddo

  end subroutine test_converge


  subroutine step(xn,xp,ye,logge,dx,dye,dxn,dxp,det_out,dxdn_out,dxdp_out,dyedn_out,dyedp_out)
    real(8),intent(in) :: xn,xp,ye
    real(8),intent(in) :: logge(n_spec)
    real(8),intent(out) :: dx,dye,dxn,dxp
    real(8),intent(out),optional :: det_out,dxdn_out,dxdp_out,dyedn_out,dyedp_out
    
    real(8) :: logx(n_spec), x(n_spec)

    real(8) :: xsum,yesum
    real(8) :: det,dxdn,dxdp,dyedn,dyedp

    real(8) :: ave_a,a11,a12,a21,a22

    logx(1:n_spec) = logge(1:n_spec) + z(1:n_spec)*xp + n(1:n_spec)*xn 

!write(6,*) logx(:)
!stop
    logx(1:n_spec) = max(-3d2,min(3d2,logx(1:n_spec)))

    x(1:n_spec) = 10d0**logx(1:n_spec)
    
    xsum = sum(x(1:n_spec))
    yesum= sum(zai(1:n_spec)*x(1:n_spec))

    dxdn = sum(n(1:n_spec)*x(1:n_spec))
    dxdp = sum(z(1:n_spec)*x(1:n_spec))
    
    dyedn = sum(n(1:n_spec)*zai(1:n_spec)*x(1:n_spec))
    dyedp = sum(z(1:n_spec)*zai(1:n_spec)*x(1:n_spec))

    dxdn = dxdn/xsum
    dxdp = dxdp/xsum
    
    dyedn = dyedn/yesum! - dxdn
    dyedp = dyedp/yesum! - dxdp

    dx  = log10(xsum)
    dye = log10(yesum/ye)

    
    det = dxdn*dyedp - dxdp*dyedn
    if( det == 0d0)then
       !write(6,*) "det = 0",dxdn*dyedp, dxdp*dyedn, dxdn,dxdp,dyedn,dyedp
       !det = 1d0
       dxn = -( dx*dyedp-dye*dxdp)
       dxp = -(-dx*dyedn+dye*dxdn)
       ! ave_a = 2d0
       ! a11 = - (x(1)+2d0*x(2))/x(3)/ave_a
       ! a12 = - (x(2)+2d0*x(1))/x(3)/ave_a
       ! a21 = - (     4d0*x(2))/x(3)/ave_a
       ! a22 = - (     2d0*x(2))/x(3)/ave_a
       ! det = ave_a**2*(a11+a22+a11*a22 - (a12+a21+a12*a21))
       
       !write(6,*) det
       !stop
    else
       dxn =-( dx*dyedp-dye*dxdp)/det
       dxp =-(-dx*dyedn+dye*dxdn)/det
    endif

    !dxn =-( dx*dyedp-dye*dxdp)/det
    !dxp =-(-dx*dyedn+dye*dxdn)/det

    if(present(det_out)) det_out = det
    if(present(dxdn_out)) dxdn_out = dxdn
    if(present(dxdp_out)) dxdp_out = dxdp
    if(present(dyedn_out)) dyedn_out = dyedn
    if(present(dyedp_out)) dyedp_out = dyedp
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

  

end module module_nse
