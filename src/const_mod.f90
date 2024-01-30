module const
!  implicit real(8) (a-h,o-z)
  implicit none
  
! ! physical constants
  real(8) :: g,h,hbar,qe,avo,clight,kerg,ev2erg,mev2erg,mev2k,kev,amu,mn,mp,me,rbohr,fine,hion!,pi

  ! parameters used for Togashi EOS
  double precision, parameter :: pi       = atan(1.d0)*4.d0
  double precision, parameter :: mevtoerg = 624150.6479963235d0
  double precision, parameter :: bm       = 1.660539274d-24      ! g
  double precision, parameter :: cc       = 29979245800.d0       ! cm/s
  double precision, parameter :: emev     = 0.51099907d0         ! MeV
  double precision, parameter :: aradmev  = 1.371594d26*mevtoerg ! MeV/cm^3/MeV^4
  double precision, parameter :: daradmev = 4.72927637d16
  double precision, parameter :: kB       = 8.6171d-5            ! eV/K
  double precision, parameter :: hp       = 4.13567d-15          ! eV s
  double precision, parameter :: hb       = hp/2.d0/pi           ! eV s
  double precision, parameter :: hc       = hp*cc*1.d-6          ! MeV cm
  double precision, parameter :: hbc      = hb*cc*1.d-6          ! MeV cm
  double precision, parameter :: bmev     = bm*cc*cc*mevtoerg    ! MeV
  double precision, parameter :: pmev     = 938.272013d0         ! MeV
  double precision, parameter :: nmev     = 939.565346d0         ! MeV
  double precision, parameter :: amev     = 931.49432d0          ! MeV
  
  real(8) :: eulercon, a2rad,rad2a
  parameter( &!pi       = 3.1415926535897932384d0, &
             eulercon = 0.577215664901532861d0, &
             a2rad    = pi/180.0d0,  rad2a = 180.0d0/pi)

  ! Togashi Parameters:
  parameter(g       = 6.6742867d-8, &
            h       = hp/mevtoerg/1.d6, &
            hbar    = 0.5d0 * h/pi, &
            qe      = 4.8032042712d-10, &
            avo     = 6.0221417930d23, &
            clight  = 2.99792458d10, &
            kerg    = kB/mevtoerg/1.d6, &
            ev2erg  = 1.60217648740d-12, &
            mev2erg = ev2erg*1.d6,      &
            mev2k   = mev2erg/kerg,     &
            kev     = kerg/ev2erg, &
            amu     = 1.66053878283d-24, &
            mn      = 1.67492721184d-24, &
            mp      = 1.67262163783d-24, &
            me      = emev/mevtoerg/cc/cc, &
            rbohr   = hbar*hbar/(me * qe * qe), &
            fine    = qe*qe/(hbar*clight), &
            hion    = 13.605698140d0)

  ! parameter(g       = 6.6742867d-8, &
  !           h       = 6.6260689633d-27, &
  !           hbar    = 0.5d0 * h/pi, &
  !           qe      = 4.8032042712d-10, &
  !           avo     = 6.0221417930d23, &
  !           clight  = 2.99792458d10, &
  !           kerg    = 1.380650424d-16, &
  !           ev2erg  = 1.60217648740d-12, &
  !           mev2erg = ev2erg*1.d6,      &
  !           mev2k   = mev2erg/kerg,     &
  !           kev     = kerg/ev2erg, &
  !           amu     = 1.66053878283d-24, &
  !           mn      = 1.67492721184d-24, &
  !           mp      = 1.67262163783d-24, &
  !           me      = 9.1093821545d-28, &
  !           rbohr   = hbar*hbar/(me * qe * qe), &
  !           fine    = qe*qe/(hbar*clight), &
  !           hion    = 13.605698140d0)
  
  real(8) :: ssol,asol,weinlam,weinfre,rhonuc
  parameter(ssol     = 2.d0*pi**5*kerg**4/(15.d0*h**3*clight**2), &
            asol     = 4.0d0 * ssol / clight, &
            weinlam  = h*clight/(kerg * 4.965114232d0), &
            weinfre  = 2.821439372d0*kerg/h, &
            rhonuc   = 2.342d14 )
  

! astronomical constants
  real(8) :: msol,rsol,lsol,mearth,rearth,ly,pc,au,secyer
  parameter(msol    = 1.9892d33, &
            rsol    = 6.95997d10, &
            lsol    = 3.8268d33, &
            mearth  = 5.9764d27, &
            rearth  = 6.37d8, &
            ly      = 9.460528d17, &
            pc      = 3.261633d0 * ly, &
            au      = 1.495978921d13, &
            secyer  = 3.1558149984d7)


  ! 
  real(8) :: na,mev,arad,cl,mu,mumev,mpmev,mnmev,mamev,memev!,kb
  parameter( na    = avo, &
             mev   = mev2erg, &
             cl    = clight, &
             mu    = amu,  &
!             kb    = kerg, &
             arad  = pi**2/(15.d0*hbar**3*cl**3)*mev2erg**4, &
             mumev = mu*cl*cl/mev2erg, &
             mpmev = mp*cl*cl/mev2erg, &
             mnmev = mn*cl*cl/mev2erg, &
             memev = me*cl*cl/mev2erg, &
             mamev = 3727.379378d0     )

      
end module const
