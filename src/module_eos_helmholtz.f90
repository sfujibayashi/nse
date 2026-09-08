module module_eos_helmholtz
  implicit none
  private
  public :: eos_get_hmin
  public :: eos_all
  public :: eos_get_eps_min
  public :: init_eos
  public :: eos_get_range
  public :: eos_get_misc
  public :: eos_get_temp_from_pres
  public :: eos_get_temp_from_eps
  public :: eos_get_eps_at_temp_min

!!! originally in const.dek
  
  ! math constants
  real*8, parameter:: pi       = 3.1415926535897932384d0, &
       eulercon = 0.577215664901532861d0, &
       a2rad    = pi/180.0d0,  rad2a = 180.0d0/pi
  
  ! physical constants
  real*8, parameter :: g       = 6.6742867d-8, &
       ! h       = 6.6260689633d-27, &
       h = 6.62607015d-27, &
       hbar    = 0.5d0 * h/pi, &
       qe      = 4.8032042712d-10, &
       avo     = 6.0221417930d23, &
       clight  = 2.99792458d10, &
       kerg    = 1.380650424d-16, &
       ! ev2erg  = 1.60217648740d-12, &
       ev2erg = 1.602176634d-12, &
       kev     = kerg/ev2erg, &
       amu     = 1.66053878283d-24, &
       mn      = 1.67492721184d-24, &
       mp      = 1.67262163783d-24, &
       me      = 9.1093821545d-28, &
       rbohr   = hbar*hbar/(me * qe * qe), &
       fine    = qe*qe/(hbar*clight), &
       hion    = 13.605698140d0, & 
       ssol    = 5.6704d-5, &
       asol    = 4.0d0 * ssol / clight, &
       weinlam = h*clight/(kerg * 4.965114232d0), &
       weinfre = 2.821439372d0*kerg/h, &
       rhonuc  = 2.342d14


  ! astronomical constants
  real*8, parameter :: msol    = 1.9892d33, &
       rsol    = 6.95997d10, &
       lsol    = 3.8268d33, &
       mearth  = 5.9764d27, &
       rearth  = 6.37d8, &
       ly      = 9.460528d17, &
       pc      = 3.261633d0 * ly, &
       au      = 1.495978921d13, &
       secyer  = 3.1558149984d7

  
!!! originally in helm_table_strage.dek
  
  ! sizes of the tables
  ! normal table, big table, bigger table, denser bigger table
  
  integer          imax,jmax

  ! original
  !      parameter        (imax = 211, jmax = 71)

  ! standard
  !      parameter        (imax = 271, jmax = 101)

  ! twice as dense
  parameter        (imax = 541, jmax = 201)

  ! half as dense
  !      parameter        (imax = 136, jmax = 51)



  ! for the electrons
  ! density and temperature
  double precision tlo,thi,tstp,tstpi,dlo,dhi,dstp,dstpi
  double precision d(imax),t(jmax)
  !      common /dttabc1/ d,t, &
  !                       tlo,thi,tstp,tstpi,dlo,dhi,dstp,dstpi

  ! for the helmholtz free energy tables
  double precision f(imax,jmax),fd(imax,jmax), &
       ft(imax,jmax),fdd(imax,jmax),ftt(imax,jmax), &
       fdt(imax,jmax),fddt(imax,jmax),fdtt(imax,jmax), &
       fddtt(imax,jmax)

  ! common /frtabc1/ f,fd, &
  !                  ft,fdd,ftt, &
  !                  fdt,fddt,fdtt, &
  !                  fddtt

  ! for the pressure derivative with density tables
  double precision dpdf(imax,jmax),dpdfd(imax,jmax), &
       dpdft(imax,jmax),dpdfdt(imax,jmax)

  ! common /dpdtab1/ dpdf,dpdfd, &
  !      dpdft,dpdfdt


  ! for chemical potential tables
  double precision ef(imax,jmax),efd(imax,jmax), &
       eft(imax,jmax),efdt(imax,jmax)

  ! common /eftabc1/ ef,efd, &
  !      eft,efdt


  ! for the number density tables
  double precision xf(imax,jmax),xfd(imax,jmax), &
       xft(imax,jmax),xfdt(imax,jmax)

  ! common /xftabc1/ xf,xfd, &
  !      xft,xfdt


  ! for storing the differences
  double precision dt_sav(jmax),dt2_sav(jmax), &
       dti_sav(jmax),dt2i_sav(jmax),dt3i_sav(jmax), &
       dd_sav(imax),dd2_sav(imax), &
       ddi_sav(imax),dd2i_sav(imax),dd3i_sav(imax)

  ! common /diftabc1/dt_sav,dt2_sav, &
  !      dti_sav,dt2i_sav,dt3i_sav, &
  !      dd_sav,dd2_sav, &
  !      ddi_sav,dd2i_sav,dd3i_sav




  ! for the ions
  ! density and temperature
  double precision tion_lo,tion_hi,tion_stp,tion_stpi, &
       dion_lo,dion_hi,dion_stp,dion_stpi
  double precision dion(imax),tion(jmax)
  ! common /dttabc2/ dion,tion, &
  !      tion_lo,tion_hi,tion_stp,tion_stpi, &
  !      dion_lo,dion_hi,dion_stp,dion_stpi

  ! for the helmholtz free energy tables
  double precision fion(imax,jmax),fiond(imax,jmax), &
       fiont(imax,jmax),fiondd(imax,jmax), &
       fiontt(imax,jmax),fiondt(imax,jmax), &
       fionddt(imax,jmax),fiondtt(imax,jmax), &
       fionddtt(imax,jmax)

  ! common /frtabc2/ fion,fiond, &
  !      fiont,fiondd,fiontt, &
  !      fiondt,fionddt,fiondtt, &
  !      fionddtt

  ! for the pressure derivative with density tables
  double precision dpiondf(imax,jmax),dpiondfd(imax,jmax), &
       dpiondft(imax,jmax),dpiondfdt(imax,jmax)

  ! common /dpdtab2/ dpiondf,dpiondfd, &
  !      dpiondft,dpiondfdt


  ! for chemical potential tables
  double precision efion(imax,jmax),efiond(imax,jmax), &
       efiont(imax,jmax),efiondt(imax,jmax)

  ! common /eftabc2/ efion,efiond, &
  !      efiont,efiondt


  ! for the number density tables
  double precision xfion(imax,jmax),xfiond(imax,jmax), &
       xfiont(imax,jmax),xfiondt(imax,jmax)

  ! common /xftabc2/ xfion,xfiond, &
  !      xfiont,xfiondt


  ! for storing the differences
  double precision dt_sav_ion(jmax),dt2_sav_ion(jmax), &
       dti_sav_ion(jmax),dt2i_sav_ion(jmax), &
       dt3i_sav_ion(jmax),dd_sav_ion(imax), &
       dd2_sav_ion(imax),ddi_sav_ion(imax), &
       dd2i_sav_ion(imax),dd3i_sav_ion(imax)

  ! common /diftabc2/dt_sav_ion,dt2_sav_ion, &
  !      dti_sav_ion,dt2i_sav_ion, &
  !      dt3i_sav_ion,dd_sav_ion, &
  !      dd2_sav_ion,ddi_sav_ion, &
  !      dd2i_sav_ion,dd3i_sav_ion



!!! originally in vector_eos.dek
  
  ! declaration for pipelining the eos routines
  
  ! maximum length of the row vector
  integer   nrowmax
  parameter (nrowmax = 1000)
  !      parameter (nrowmax = 10000)


  ! maximum number of isotopes
  integer   irowmax
  !      parameter (irowmax = 1)
  parameter (irowmax = 100)


  ! maximum number of ionization stages
  integer   jstagemax
  !      parameter (jstagemax = 1)
  parameter (jstagemax = 30)


  ! failure of an eos
  logical ::  eosfail = .false.
  ! common /eosfc1/  eosfail



  ! lower and upper limits of the loop over rows
  integer          jlo_eos,jhi_eos
  ! common /eosvec2/ jlo_eos,jhi_eos



  ! thermodynamic and composition inputs
  double precision &
       temp_row(nrowmax),den_row(nrowmax), &
       abar_row(nrowmax),zbar_row(nrowmax), &
       zeff_row(nrowmax),ye_row(nrowmax)

  ! common /thinp/ &
  !           temp_row,den_row, &
  !           abar_row,zbar_row, &
  !           zeff_row,ye_row


  ! composition input
  integer          niso
  double precision xmass_row(irowmax,nrowmax), &
       aion_row(irowmax,nrowmax), &
       zion_row(irowmax,nrowmax)
  ! common  /cmpinp/ xmass_row,aion_row,zion_row,niso



  ! composition output
  double precision frac_row(jstagemax,irowmax,nrowmax)
  ! common /cmpout/  frac_row


  ! composition output for sneos
  double precision xn_row(nrowmax),xp_row(nrowmax), &
       xa_row(nrowmax),xhv_row(nrowmax), &
       xmuhat_row(nrowmax)
  ! common /cmpout2/ xn_row,xp_row, &
  !      xa_row,xhv_row, &
  !      xmuhat_row

  


  ! totals and their derivatives
  double precision &
       ptot_row(nrowmax), &
       dpt_row(nrowmax),dpd_row(nrowmax), &
       dpa_row(nrowmax),dpz_row(nrowmax), &
       dpdd_row(nrowmax),dpdt_row(nrowmax), &
       dpda_row(nrowmax),dpdz_row(nrowmax), &
       dptt_row(nrowmax),dpta_row(nrowmax), &
       dptz_row(nrowmax),dpaa_row(nrowmax), &
       dpaz_row(nrowmax),dpzz_row(nrowmax)

  ! common /ptotc1/ &
  !      ptot_row, &
  !      dpt_row,dpd_row, &
  !      dpa_row,dpz_row, &
  !      dpdd_row,dpdt_row, &
  !      dpda_row,dpdz_row, &
  !      dptt_row,dpta_row, &
  !      dptz_row,dpaa_row, &
  !      dpaz_row,dpzz_row

  double precision &
       etot_row(nrowmax), &
       det_row(nrowmax),ded_row(nrowmax), &
       dea_row(nrowmax),dez_row(nrowmax), &
       dedd_row(nrowmax),dedt_row(nrowmax), &
       deda_row(nrowmax),dedz_row(nrowmax), &
       dett_row(nrowmax),deta_row(nrowmax), &
       detz_row(nrowmax),deaa_row(nrowmax), &
       deaz_row(nrowmax),dezz_row(nrowmax)

  ! common /etotc1/ &
  !      etot_row, &
  !      det_row,ded_row, &
  !      dea_row,dez_row, &
  !      dedd_row,dedt_row, &
  !      deda_row,dedz_row, &
  !      dett_row,deta_row, &
  !      detz_row,deaa_row, &
  !      deaz_row,dezz_row

  double precision &
       stot_row(nrowmax), &
       dst_row(nrowmax),dsd_row(nrowmax), &
       dsa_row(nrowmax),dsz_row(nrowmax), &
       dsdd_row(nrowmax),dsdt_row(nrowmax), &
       dsda_row(nrowmax),dsdz_row(nrowmax), &
       dstt_row(nrowmax),dsta_row(nrowmax), &
       dstz_row(nrowmax),dsaa_row(nrowmax), &
       dsaz_row(nrowmax),dszz_row(nrowmax)

  ! common /stotc1/ &
  !      stot_row, &
  !      dst_row,dsd_row, &
  !      dsa_row,dsz_row, &
  !      dsdd_row,dsdt_row, &
  !      dsda_row,dsdz_row, &
  !      dstt_row,dsta_row, &
  !      dstz_row,dsaa_row, &
  !      dsaz_row,dszz_row



  ! radiation contributions
  double precision &
       prad_row(nrowmax), &
       dpradt_row(nrowmax),dpradd_row(nrowmax), &
       dprada_row(nrowmax),dpradz_row(nrowmax), &
       dpraddd_row(nrowmax),dpraddt_row(nrowmax), &
       dpradda_row(nrowmax),dpraddz_row(nrowmax), &
       dpradtt_row(nrowmax),dpradta_row(nrowmax), &
       dpradtz_row(nrowmax),dpradaa_row(nrowmax), &
       dpradaz_row(nrowmax),dpradzz_row(nrowmax)
  ! common /thprad/ &
  !      prad_row, &
  !      dpradt_row,dpradd_row, &
  !      dprada_row,dpradz_row, &
  !      dpraddd_row,dpraddt_row, &
  !      dpradda_row,dpraddz_row, &
  !      dpradtt_row,dpradta_row, &
  !      dpradtz_row,dpradaa_row, &
  !      dpradaz_row,dpradzz_row


  double precision &
       erad_row(nrowmax), &
       deradt_row(nrowmax),deradd_row(nrowmax), &
       derada_row(nrowmax),deradz_row(nrowmax), &
       deraddd_row(nrowmax),deraddt_row(nrowmax), &
       deradda_row(nrowmax),deraddz_row(nrowmax), &
       deradtt_row(nrowmax),deradta_row(nrowmax), &
       deradtz_row(nrowmax),deradaa_row(nrowmax), &
       deradaz_row(nrowmax),deradzz_row(nrowmax)
  ! common /therad/ &
  !      erad_row, &
  !      deradt_row,deradd_row, &
  !      derada_row,deradz_row, &
  !      deraddd_row,deraddt_row, &
  !      deradda_row,deraddz_row, &
  !      deradtt_row,deradta_row, &
  !      deradtz_row,deradaa_row, &
  !      deradaz_row,deradzz_row


  double precision &
       srad_row(nrowmax), &
       dsradt_row(nrowmax),dsradd_row(nrowmax), &
       dsrada_row(nrowmax),dsradz_row(nrowmax), &
       dsraddd_row(nrowmax),dsraddt_row(nrowmax), &
       dsradda_row(nrowmax),dsraddz_row(nrowmax), &
       dsradtt_row(nrowmax),dsradta_row(nrowmax), &
       dsradtz_row(nrowmax),dsradaa_row(nrowmax), &
       dsradaz_row(nrowmax),dsradzz_row(nrowmax)
  ! common /thsrad/ &
  !      srad_row, &
  !      dsradt_row,dsradd_row, &
  !      dsrada_row,dsradz_row, &
  !      dsraddd_row,dsraddt_row, &
  !      dsradda_row,dsraddz_row, &
  !      dsradtt_row,dsradta_row, &
  !      dsradtz_row,dsradaa_row, &
  !      dsradaz_row,dsradzz_row



  ! gas contributions
  double precision &
       pgas_row(nrowmax), &
       dpgast_row(nrowmax),dpgasd_row(nrowmax), &
       dpgasa_row(nrowmax),dpgasz_row(nrowmax), &
       dpgasdd_row(nrowmax),dpgasdt_row(nrowmax), &
       dpgasda_row(nrowmax),dpgasdz_row(nrowmax), &
       dpgastt_row(nrowmax),dpgasta_row(nrowmax), &
       dpgastz_row(nrowmax),dpgasaa_row(nrowmax), &
       dpgasaz_row(nrowmax),dpgaszz_row(nrowmax)

  ! common /thpgasc1/ &
  !      pgas_row, &
  !      dpgast_row,dpgasd_row, &
  !      dpgasa_row,dpgasz_row, &
  !      dpgasdd_row,dpgasdt_row, &
  !      dpgasda_row,dpgasdz_row, &
  !      dpgastt_row,dpgasta_row, &
  !      dpgastz_row,dpgasaa_row, &
  !      dpgasaz_row,dpgaszz_row

  double precision &
       egas_row(nrowmax), &
       degast_row(nrowmax),degasd_row(nrowmax), &
       degasa_row(nrowmax),degasz_row(nrowmax), &
       degasdd_row(nrowmax),degasdt_row(nrowmax), &
       degasda_row(nrowmax),degasdz_row(nrowmax), &
       degastt_row(nrowmax),degasta_row(nrowmax), &
       degastz_row(nrowmax),degasaa_row(nrowmax), &
       degasaz_row(nrowmax),degaszz_row(nrowmax)

  ! common /thegasc2/ &
  !      egas_row, &
  !      degast_row,degasd_row, &
  !      degasa_row,degasz_row, &
  !      degasdd_row,degasdt_row, &
  !      degasda_row,degasdz_row, &
  !      degastt_row,degasta_row, &
  !      degastz_row,degasaa_row, &
  !      degasaz_row,degaszz_row

  double precision &
       sgas_row(nrowmax), &
       dsgast_row(nrowmax),dsgasd_row(nrowmax), &
       dsgasa_row(nrowmax),dsgasz_row(nrowmax), &
       dsgasdd_row(nrowmax),dsgasdt_row(nrowmax), &
       dsgasda_row(nrowmax),dsgasdz_row(nrowmax), &
       dsgastt_row(nrowmax),dsgasta_row(nrowmax), &
       dsgastz_row(nrowmax),dsgasaa_row(nrowmax), &
       dsgasaz_row(nrowmax),dsgaszz_row(nrowmax)

  ! common /thsgasc1/ &
  !      sgas_row, &
  !      dsgast_row,dsgasd_row, &
  !      dsgasa_row,dsgasz_row, &
  !      dsgasdd_row,dsgasdt_row, &
  !      dsgasda_row,dsgasdz_row, &
  !      dsgastt_row,dsgasta_row, &
  !      dsgastz_row,dsgasaa_row, &
  !      dsgasaz_row,dsgaszz_row






  ! ion contributions
  double precision &
       pion_row(nrowmax), &
       dpiont_row(nrowmax),dpiond_row(nrowmax), &
       dpiona_row(nrowmax),dpionz_row(nrowmax), &
       dpiondd_row(nrowmax),dpiondt_row(nrowmax), &
       dpionda_row(nrowmax),dpiondz_row(nrowmax), &
       dpiontt_row(nrowmax),dpionta_row(nrowmax), &
       dpiontz_row(nrowmax),dpionaa_row(nrowmax), &
       dpionaz_row(nrowmax),dpionzz_row(nrowmax)
  ! common /thpion/ &
  !      pion_row, &
  !      dpiont_row,dpiond_row, &
  !      dpiona_row,dpionz_row, &
  !      dpiondd_row,dpiondt_row, &
  !      dpionda_row,dpiondz_row, &
  !      dpiontt_row,dpionta_row, &
  !      dpiontz_row,dpionaa_row, &
  !      dpionaz_row,dpionzz_row


  double precision &
       eion_row(nrowmax), &
       deiont_row(nrowmax),deiond_row(nrowmax), &
       deiona_row(nrowmax),deionz_row(nrowmax), &
       deiondd_row(nrowmax),deiondt_row(nrowmax), &
       deionda_row(nrowmax),deiondz_row(nrowmax), &
       deiontt_row(nrowmax),deionta_row(nrowmax), &
       deiontz_row(nrowmax),deionaa_row(nrowmax), &
       deionaz_row(nrowmax),deionzz_row(nrowmax)
  ! common /theion/ &
  !      eion_row, &
  !      deiont_row,deiond_row, &
  !      deiona_row,deionz_row, &
  !      deiondd_row,deiondt_row, &
  !      deionda_row,deiondz_row, &
  !      deiontt_row,deionta_row, &
  !      deiontz_row,deionaa_row, &
  !      deionaz_row,deionzz_row


  double precision &
       sion_row(nrowmax), &
       dsiont_row(nrowmax),dsiond_row(nrowmax), &
       dsiona_row(nrowmax),dsionz_row(nrowmax), &
       dsiondd_row(nrowmax),dsiondt_row(nrowmax), &
       dsionda_row(nrowmax),dsiondz_row(nrowmax), &
       dsiontt_row(nrowmax),dsionta_row(nrowmax), &
       dsiontz_row(nrowmax),dsionaa_row(nrowmax), &
       dsionaz_row(nrowmax),dsionzz_row(nrowmax)
  ! common /thsion/ &
  !      sion_row, &
  !      dsiont_row,dsiond_row, &
  !      dsiona_row,dsionz_row, &
  !      dsiondd_row,dsiondt_row, &
  !      dsionda_row,dsiondz_row, &
  !      dsiontt_row,dsionta_row, &
  !      dsiontz_row,dsionaa_row, &
  !      dsionaz_row,dsionzz_row

  double precision &
       etaion_row(nrowmax), &
       detait_row(nrowmax),detaid_row(nrowmax), &
       detaia_row(nrowmax),detaiz_row(nrowmax), &
       detaidd_row(nrowmax),detaidt_row(nrowmax), &
       detaida_row(nrowmax),detaidz_row(nrowmax), &
       detaitt_row(nrowmax),detaita_row(nrowmax), &
       detaitz_row(nrowmax),detaiaa_row(nrowmax), &
       detaiaz_row(nrowmax),detaizz_row(nrowmax)
  ! common /thetaion/ &
  !      etaion_row, &
  !      detait_row,detaid_row, &
  !      detaia_row,detaiz_row, &
  !      detaidd_row,detaidt_row, &
  !      detaida_row,detaidz_row, &
  !      detaitt_row,detaita_row, &
  !      detaitz_row,detaiaa_row, &
  !      detaiaz_row,detaizz_row


  double precision &
       xni_row(nrowmax),xnim_row(nrowmax), &
       dxnit_row(nrowmax),dxnid_row(nrowmax), &
       dxnia_row(nrowmax),dxniz_row(nrowmax), &
       dxnidd_row(nrowmax),dxnidt_row(nrowmax), &
       dxnida_row(nrowmax),dxnidz_row(nrowmax), &
       dxnitt_row(nrowmax),dxnita_row(nrowmax), &
       dxnitz_row(nrowmax),dxniaa_row(nrowmax), &
       dxniaz_row(nrowmax),dxnizz_row(nrowmax)
  ! common /th_xni_ion/ &
  !      xni_row,xnim_row, &
  !      dxnit_row,dxnid_row, &
  !      dxnia_row,dxniz_row, &
  !      dxnidd_row,dxnidt_row, &
  !      dxnida_row,dxnidz_row, &
  !      dxnitt_row,dxnita_row, &
  !      dxnitz_row,dxniaa_row, &
  !      dxniaz_row,dxnizz_row



  ! electron-positron contributions

  double precision &
       etaele_row(nrowmax),etapos_row(nrowmax), &
       detat_row(nrowmax),detad_row(nrowmax), &
       detaa_row(nrowmax),detaz_row(nrowmax), &
       detadd_row(nrowmax),detadt_row(nrowmax), &
       detada_row(nrowmax),detadz_row(nrowmax), &
       detatt_row(nrowmax),detata_row(nrowmax), &
       detatz_row(nrowmax),detaaa_row(nrowmax), &
       detaaz_row(nrowmax),detazz_row(nrowmax)

  ! common /etapc1/ &
  !      etaele_row,etapos_row, &
  !      detat_row,detad_row, &
  !      detaa_row,detaz_row, &
  !      detadd_row,detadt_row, &
  !      detada_row,detadz_row, &
  !      detatt_row,detata_row, &
  !      detatz_row,detaaa_row, &
  !      detaaz_row,detazz_row

  double precision &
       pele_row(nrowmax),ppos_row(nrowmax), &
       dpept_row(nrowmax),dpepd_row(nrowmax), &
       dpepa_row(nrowmax),dpepz_row(nrowmax), &
       dpepdd_row(nrowmax),dpepdt_row(nrowmax), &
       dpepda_row(nrowmax),dpepdz_row(nrowmax), &
       dpeptt_row(nrowmax),dpepta_row(nrowmax), &
       dpeptz_row(nrowmax),dpepaa_row(nrowmax), &
       dpepaz_row(nrowmax),dpepzz_row(nrowmax)

  ! common /thpepc1/ &
  !      pele_row,ppos_row, &
  !      dpept_row,dpepd_row, &
  !      dpepa_row,dpepz_row, &
  !      dpepdd_row,dpepdt_row, &
  !      dpepda_row,dpepdz_row, &
  !      dpeptt_row,dpepta_row, &
  !      dpeptz_row,dpepaa_row, &
  !      dpepaz_row,dpepzz_row


  double precision &
       eele_row(nrowmax),epos_row(nrowmax), &
       deept_row(nrowmax),deepd_row(nrowmax), &
       deepa_row(nrowmax),deepz_row(nrowmax), &
       deepdd_row(nrowmax),deepdt_row(nrowmax), &
       deepda_row(nrowmax),deepdz_row(nrowmax), &
       deeptt_row(nrowmax),deepta_row(nrowmax), &
       deeptz_row(nrowmax),deepaa_row(nrowmax), &
       deepaz_row(nrowmax),deepzz_row(nrowmax)

  ! common /theepc1/ &
  !      eele_row,epos_row, &
  !      deept_row,deepd_row, &
  !      deepa_row,deepz_row, &
  !      deepdd_row,deepdt_row, &
  !      deepda_row,deepdz_row, &
  !      deeptt_row,deepta_row, &
  !      deeptz_row,deepaa_row, &
  !      deepaz_row,deepzz_row


  double precision &
       sele_row(nrowmax),spos_row(nrowmax), &
       dsept_row(nrowmax),dsepd_row(nrowmax), &
       dsepa_row(nrowmax),dsepz_row(nrowmax), &
       dsepdd_row(nrowmax),dsepdt_row(nrowmax), &
       dsepda_row(nrowmax),dsepdz_row(nrowmax), &
       dseptt_row(nrowmax),dsepta_row(nrowmax), &
       dseptz_row(nrowmax),dsepaa_row(nrowmax), &
       dsepaz_row(nrowmax),dsepzz_row(nrowmax)

  ! common /thsepc1/ &
  !      sele_row,spos_row, &
  !      dsept_row,dsepd_row, &
  !      dsepa_row,dsepz_row, &
  !      dsepdd_row,dsepdt_row, &
  !      dsepda_row,dsepdz_row, &
  !      dseptt_row,dsepta_row, &
  !      dseptz_row,dsepaa_row, &
  !      dsepaz_row,dsepzz_row


  double precision &
       xne_row(nrowmax),xnp_row(nrowmax),xnem_row(nrowmax), &
       dxnet_row(nrowmax),dxned_row(nrowmax), &
       dxnea_row(nrowmax),dxnez_row(nrowmax), &
       dxnedd_row(nrowmax),dxnedt_row(nrowmax), &
       dxneda_row(nrowmax),dxnedz_row(nrowmax), &
       dxnett_row(nrowmax),dxneta_row(nrowmax), &
       dxnetz_row(nrowmax),dxneaa_row(nrowmax), &
       dxneaz_row(nrowmax),dxnezz_row(nrowmax)

  ! common /thxnec1/ &
  !      xne_row,xnp_row,xnem_row, &
  !      dxnet_row,dxned_row, &
  !      dxnea_row,dxnez_row, &
  !      dxnedd_row,dxnedt_row, &
  !      dxneda_row,dxnedz_row, &
  !      dxnett_row,dxneta_row, &
  !      dxnetz_row,dxneaa_row, &
  !      dxneaz_row,dxnezz_row



  ! ionization potential contributions
  double precision pip_row(nrowmax), &
       dpipt_row(nrowmax), dpipd_row(nrowmax), &
       dpipa_row(nrowmax), dpipz_row(nrowmax), &
       eip_row(nrowmax), & 
       deipt_row(nrowmax), deipd_row(nrowmax), &
       deipa_row(nrowmax), deipz_row(nrowmax), &
       sip_row(nrowmax), &
       dsipt_row(nrowmax), dsipd_row(nrowmax), &
       dsipa_row(nrowmax), dsipz_row(nrowmax)
  ! common /thxip/   pip_row,dpipt_row,dpipd_row,dpipa_row,dpipz_row, &
  !      eip_row,deipt_row,deipd_row,deipa_row,deipz_row, &
  !      sip_row,dsipt_row,dsipd_row,dsipa_row,dsipz_row



  ! coulomb contributions
  double precision &
       pcou_row(nrowmax), &
       dpcout_row(nrowmax),dpcoud_row(nrowmax), &
       dpcoua_row(nrowmax),dpcouz_row(nrowmax), &
       ecou_row(nrowmax), &
       decout_row(nrowmax),decoud_row(nrowmax), &
       decoua_row(nrowmax),decouz_row(nrowmax), &
       scou_row(nrowmax), &
       dscout_row(nrowmax),dscoud_row(nrowmax), &
       dscoua_row(nrowmax),dscouz_row(nrowmax), &
       plasg_row(nrowmax)
  ! common /thcou/ &
  !      pcou_row, &
  !      dpcout_row,dpcoud_row, &
  !      dpcoua_row,dpcouz_row, &
  !      ecou_row, &
  !      decout_row,decoud_row, &
  !      decoua_row,decouz_row, &
  !      scou_row, &
  !      dscout_row,dscoud_row, &
  !      dscoua_row,dscouz_row, &
  !      plasg_row


  ! thermodynamic consistency checks; maxwell relations
  double precision &
       dse_row(nrowmax),dpe_row(nrowmax),dsp_row(nrowmax)
  ! common /thmax/ &
  !      dse_row,dpe_row,dsp_row


  ! derivative based quantities for the gas
  double precision &
       cp_gas_row(nrowmax), &
       dcp_gasdd_row(nrowmax),dcp_gasdt_row(nrowmax), &
       dcp_gasda_row(nrowmax),dcp_gasdz_row(nrowmax), &
       cv_gas_row(nrowmax), &
       dcv_gasdd_row(nrowmax),dcv_gasdt_row(nrowmax), &
       dcv_gasda_row(nrowmax),dcv_gasdz_row(nrowmax)

  ! common /thdergc1/ &
  !      cp_gas_row, &
  !      dcp_gasdd_row,dcp_gasdt_row, &
  !      dcp_gasda_row,dcp_gasdz_row, &
  !      cv_gas_row, &
  !      dcv_gasdd_row,dcv_gasdt_row, &
  !      dcv_gasda_row,dcv_gasdz_row

  double precision &
       gam1_gas_row(nrowmax), &
       dgam1_gasdd_row(nrowmax),dgam1_gasdt_row(nrowmax), &
       dgam1_gasda_row(nrowmax),dgam1_gasdz_row(nrowmax), &
       gam2_gas_row(nrowmax), &
       dgam2_gasdd_row(nrowmax),dgam2_gasdt_row(nrowmax), &
       dgam2_gasda_row(nrowmax),dgam2_gasdz_row(nrowmax), &
       gam3_gas_row(nrowmax), &
       dgam3_gasdd_row(nrowmax),dgam3_gasdt_row(nrowmax), &
       dgam3_gasda_row(nrowmax),dgam3_gasdz_row(nrowmax), &
       nabad_gas_row(nrowmax), &
       dnab_gasdd_row(nrowmax),dnab_gasdt_row(nrowmax), &
       dnab_gasda_row(nrowmax),dnab_gasdz_row(nrowmax), &
       cs_gas_row(nrowmax), &
       dcs_gasdd_row(nrowmax),dcs_gasdt_row(nrowmax), &
       dcs_gasda_row(nrowmax),dcs_gasdz_row(nrowmax)

  ! common /thdergc2/ &
  !      gam1_gas_row, &
  !      dgam1_gasdd_row,dgam1_gasdt_row, &
  !      dgam1_gasda_row,dgam1_gasdz_row, &
  !      gam2_gas_row, &
  !      dgam2_gasdd_row,dgam2_gasdt_row, &
  !      dgam2_gasda_row,dgam2_gasdz_row, &
  !      gam3_gas_row, &
  !      dgam3_gasdd_row,dgam3_gasdt_row, &
  !      dgam3_gasda_row,dgam3_gasdz_row, &
  !      nabad_gas_row, &
  !      dnab_gasdd_row,dnab_gasdt_row, &
  !      dnab_gasda_row,dnab_gasdz_row, &
  !      cs_gas_row, &
  !      dcs_gasdd_row,dcs_gasdt_row, &
  !      dcs_gasda_row,dcs_gasdz_row



  ! derivative based quantities for the totals
  double precision &
       cp_row(nrowmax), &
       dcpdd_row(nrowmax),dcpdt_row(nrowmax), &
       dcpda_row(nrowmax),dcpdz_row(nrowmax), &
       cv_row(nrowmax), &
       dcvdd_row(nrowmax),dcvdt_row(nrowmax), &
       dcvda_row(nrowmax),dcvdz_row(nrowmax)

  ! common /thdertc1/ &
  !      cp_row, &
  !      dcpdd_row,dcpdt_row, &
  !      dcpda_row,dcpdz_row, &
  !      cv_row, &
  !      dcvdd_row,dcvdt_row, &
  !      dcvda_row,dcvdz_row

  double precision &
       gam1_row(nrowmax), &
       dgam1dd_row(nrowmax),dgam1dt_row(nrowmax), &
       dgam1da_row(nrowmax),dgam1dz_row(nrowmax), &
       gam2_row(nrowmax), &
       dgam2dd_row(nrowmax),dgam2dt_row(nrowmax), &
       dgam2da_row(nrowmax),dgam2dz_row(nrowmax), &
       gam3_row(nrowmax), &
       dgam3dd_row(nrowmax),dgam3dt_row(nrowmax), &
       dgam3da_row(nrowmax),dgam3dz_row(nrowmax), &
       nabad_row(nrowmax), &
       dnabdd_row(nrowmax),dnabdt_row(nrowmax), &
       dnabda_row(nrowmax),dnabdz_row(nrowmax), &
       cs_row(nrowmax), &
       dcsdd_row(nrowmax),dcsdt_row(nrowmax), &
       dcsda_row(nrowmax),dcsdz_row(nrowmax)

  ! common /thdertc2/ &
  !      gam1_row, &
  !      dgam1dd_row,dgam1dt_row, &
  !      dgam1da_row,dgam1dz_row, &
  !      gam2_row, &
  !      dgam2dd_row,dgam2dt_row, &
  !      dgam2da_row,dgam2dz_row, &
  !      gam3_row, &
  !      dgam3dd_row,dgam3dt_row, &
  !      dgam3da_row,dgam3dz_row, &
  !      nabad_row, &
  !      dnabdd_row,dnabdt_row, &
  !      dnabda_row,dnabdz_row, &
  !      cs_row, &
  !      dcsdd_row,dcsdt_row, &
  !      dcsda_row,dcsdz_row




  ! a few work arrays
  double precision eoswrk01(nrowmax),eoswrk02(nrowmax), &
       eoswrk03(nrowmax),eoswrk04(nrowmax)
  ! common /deedoo/  eoswrk01,eoswrk02,eoswrk03,eoswrk04


  ! for debugging
  double precision &
       crp_row(nrowmax), &
       dcrpt_row(nrowmax),dcrpd_row(nrowmax), &
       dcrpa_row(nrowmax),dcrpz_row(nrowmax), &
       dcrpdd_row(nrowmax),dcrpdt_row(nrowmax), &
       dcrpda_row(nrowmax),dcrpdz_row(nrowmax), &
       dcrptt_row(nrowmax),dcrpta_row(nrowmax), &
       dcrptz_row(nrowmax),dcrpaa_row(nrowmax), &
       dcrpaz_row(nrowmax),dcrpzz_row(nrowmax)


  ! common /crpc1/ &
  !      crp_row, &
  !      dcrpt_row,dcrpd_row, &
  !      dcrpa_row,dcrpz_row, &
  !      dcrpdd_row,dcrpdt_row, &
  !      dcrpda_row,dcrpdz_row, &
  !      dcrptt_row,dcrpta_row, &
  !      dcrptz_row,dcrpaa_row, &
  !      dcrpaz_row,dcrpzz_row




  !!! additional variables
  character(256) :: fn_eos

  real(8),parameter, public :: m_u_mev = amu*clight**2 / ev2erg / 1d6
  real(8),parameter :: m_u_g = amu
  real(8),parameter :: m_e_g = me
  real(8),parameter :: m_e_erg = me*clight**2

  real(8),parameter :: lambda_e = hbar/(m_e_g*clight), lambda_e3 = lambda_e**3
  real(8),parameter :: mev_to_erg = ev2erg*1d6

  real(8) :: temp_eos_min, temp_eos_max
contains

  subroutine init_eos(fn_eos_in)
    character(*),intent(in)  :: fn_eos_in

    fn_eos = fn_eos_in
    call read_helm_table

    temp_eos_min = t(1)
    temp_eos_max = t(jmax)
  end subroutine init_eos


  subroutine eos_all(rho_in, temp_in, ye_in, ytot_in, mexc_in, &
       eps_out, pres_out, cs2_out, entr_out, chi_out, kappa_out)
    real(8),intent(in) :: rho_in, temp_in, ye_in, ytot_in, mexc_in
    real(8),intent(out) :: eps_out, pres_out, cs2_out, entr_out
    real(8),intent(out),optional :: chi_out, kappa_out

    ! declare
    integer          i,j
    double precision temp,den,abar,zbar,ytot1,ye, &
         x,y,zz,zzi,deni,tempi,xni,dxnidd,dxnida, &
         dpepdt,dpepdd,deepdt,deepdd,dsepdd,dsepdt, &
         dpraddd,dpraddt,deraddd,deraddt,dpiondd,dpiondt, &
         deiondd,deiondt,dsraddd,dsraddt,dsiondd,dsiondt, &
         dse,dpe,dsp,kt,ktinv,prad,erad,srad,pion,eion, &
         sion,xnem,pele,eele,sele,pres,ener,entr,dpresdd, &
         dpresdt,denerdd,denerdt,dentrdd,dentrdt,cv,cp, &
         gam1,gam2,gam3,chit,chid,nabad,sound,etaele, &
         detadt,detadd,xnefer,dxnedt,dxnedd,s

    double precision pgas,dpgasdd,dpgasdt,dpgasda,dpgasdz, &
         egas,degasdd,degasdt,degasda,degasdz, &
         sgas,dsgasdd,dsgasdt,dsgasda,dsgasdz, &
         cv_gas,cp_gas,gam1_gas,gam2_gas,gam3_gas, &
         chit_gas,chid_gas,nabad_gas,sound_gas


    double precision sioncon,forth,forpi,kergavo,ikavo,asoli3,light2
    parameter        (sioncon = (2.0d0 * pi * amu * kerg)/(h*h), &
         forth   = 4.0d0/3.0d0, &
         forpi   = 4.0d0 * pi, &
         kergavo = kerg * avo, &
         ikavo   = 1.0d0/kergavo, &
         asoli3  = asol/3.0d0, &
         light2  = clight * clight)

    ! for the abar derivatives
    double precision dpradda,deradda,dsradda, &
         dpionda,deionda,dsionda, &
         dpepda,deepda,dsepda, &
         dpresda,denerda,dentrda, &
         detada,dxneda

    ! for the zbar derivatives
    double precision dpraddz,deraddz,dsraddz, &
         dpiondz,deiondz,dsiondz, &
         dpepdz,deepdz,dsepdz, &
         dpresdz,denerdz,dentrdz, &
         detadz,dxnedz

    ! for the interpolations
    integer          iat,jat
    double precision free,df_d,df_t,df_dd,df_tt,df_dt
    double precision xt,xd,mxt,mxd, &
         si0t,si1t,si2t,si0mt,si1mt,si2mt, &
         si0d,si1d,si2d,si0md,si1md,si2md, &
         dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt, &
         dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md, &
         ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt, &
         ddsi0d,ddsi1d,ddsi2d,ddsi0md,ddsi1md,ddsi2md, &
         z,psi0,dpsi0,ddpsi0,psi1,dpsi1,ddpsi1,psi2, &
         dpsi2,ddpsi2,din,h5,fi(36), &
         xpsi0,xdpsi0,xpsi1,xdpsi1,h3, &
         w0t,w1t,w2t,w0mt,w1mt,w2mt, &
         w0d,w1d,w2d,w0md,w1md,w2md


    ! for the uniform background coulomb correction
    double precision dsdd,dsda,lami,inv_lami,lamida,lamidd, &
         plasg,plasgdd,plasgdt,plasgda,plasgdz, &
         ecoul,decouldd,decouldt,decoulda,decouldz, &
         pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz, &
         scoul,dscouldd,dscouldt,dscoulda,dscouldz, &
         a1,b1,c1,d1,e1,a2,b2,c2,third,esqu
    parameter        (a1    = -0.898004d0, &
         b1    =  0.96786d0, &
         c1    =  0.220703d0, &
         d1    = -0.86097d0, &
         e1    =  2.5269d0, &
         a2    =  0.29561d0, &
         b2    =  1.9885d0, &
         c2    =  0.288675d0, &
         third =  1.0d0/3.0d0, &
         esqu  =  qe * qe)

    real(8) :: eps
    
    ! quintic hermite polynomial statement functions
    ! psi0 and its derivatives
    psi0(z)   = z**3 * ( z * (-6.0d0*z + 15.0d0) -10.0d0) + 1.0d0
    dpsi0(z)  = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
    ddpsi0(z) = z* ( z*( -120.0d0*z + 180.0d0) -60.0d0)


    ! psi1 and its derivatives
    psi1(z)   = z* ( z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
    dpsi1(z)  = z*z * ( z * (-15.0d0*z + 32.0d0) - 18.0d0) +1.0d0
    ddpsi1(z) = z * (z * (-60.0d0*z + 96.0d0) -36.0d0)


    ! psi2  and its derivatives
    psi2(z)   = 0.5d0*z*z*( z* ( z * (-z + 3.0d0) - 3.0d0) + 1.0d0)
    dpsi2(z)  = 0.5d0*z*( z*(z*(-5.0d0*z + 12.0d0) - 9.0d0) + 2.0d0)
    ddpsi2(z) = 0.5d0*(z*( z * (-20.0d0*z + 36.0d0) - 18.0d0) + 2.0d0)


    ! biquintic hermite polynomial statement function
    h5(i,j,w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md)= &
         fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t &
         + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt &
         + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t &
         + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt &
         + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t &
         + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt &
         + fi(13) *w1d*w0t   + fi(14) *w1md*w0t &
         + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt &
         + fi(17) *w2d*w0t   + fi(18) *w2md*w0t &
         + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt &
         + fi(21) *w1d*w1t   + fi(22) *w1md*w1t &
         + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt &
         + fi(25) *w2d*w1t   + fi(26) *w2md*w1t &
         + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt &
         + fi(29) *w1d*w2t   + fi(30) *w1md*w2t &
         + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt &
         + fi(33) *w2d*w2t   + fi(34) *w2md*w2t &
         + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt



    ! cubic hermite polynomial statement functions
    ! psi0 & derivatives
    xpsi0(z)  = z * z * (2.0d0*z - 3.0d0) + 1.0
    xdpsi0(z) = z * (6.0d0*z - 6.0d0)


    ! psi1 & derivatives
    xpsi1(z)  = z * ( z * (z - 2.0d0) + 1.0d0)
    xdpsi1(z) = z * (3.0d0*z - 4.0d0) + 1.0d0


    ! bicubic hermite polynomial statement function
    h3(i,j,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) = &
         fi(1)  *w0d*w0t   +  fi(2)  *w0md*w0t &
         + fi(3)  *w0d*w0mt  +  fi(4)  *w0md*w0mt &
         + fi(5)  *w0d*w1t   +  fi(6)  *w0md*w1t &
         + fi(7)  *w0d*w1mt  +  fi(8)  *w0md*w1mt &
         + fi(9)  *w1d*w0t   +  fi(10) *w1md*w0t &
         + fi(11) *w1d*w0mt  +  fi(12) *w1md*w0mt &
         + fi(13) *w1d*w1t   +  fi(14) *w1md*w1t &
         + fi(15) *w1d*w1mt  +  fi(16) *w1md*w1mt



    ! popular format statements
01  format(1x,5(a,1pe11.3))
02  format(1x,a,1p4e16.8)
03  format(1x,4(a,1pe11.3))
04  format(1x,4(a,i4))

    ! start of pipeline loop, normal execution starts here
    ! eosfail = .false.

    !       if (temp_row(j) .le. 0.0) stop 'temp less than 0 in helmeos'
    !       if (den_row(j)  .le. 0.0) stop 'den less than 0 in helmeos'

    temp  = temp_in
    den   = rho_in
    abar  = 1d0/ytot_in
    zbar  = abar*ye_in
    ytot1 = 1.0d0/abar
    ye    = ye_in

    ! initialize
    deni    = 1.0d0/den
    tempi   = 1.0d0/temp
    kt      = kerg * temp
    ktinv   = 1.0d0/kt


    ! radiation section:
    prad    = asoli3 * temp * temp * temp * temp
    dpraddd = 0.0d0
    dpraddt = 4.0d0 * prad*tempi
    dpradda = 0.0d0
    dpraddz = 0.0d0

    erad    = 3.0d0 * prad*deni
    deraddd = -erad*deni
    deraddt = 3.0d0 * dpraddt*deni
    deradda = 0.0d0
    deraddz = 0.0d0

    srad    = (prad*deni + erad)*tempi
    dsraddd = (dpraddd*deni - prad*deni*deni + deraddd)*tempi
    dsraddt = (dpraddt*deni + deraddt - srad)*tempi
    dsradda = 0.0d0
    dsraddz = 0.0d0


    ! ion section:
    xni     = avo * ytot1 * den
    dxnidd  = avo * ytot1
    dxnida  = -xni * ytot1

    pion    = xni * kt
    dpiondd = dxnidd * kt
    dpiondt = xni * kerg
    dpionda = dxnida * kt
    dpiondz = 0.0d0

    eion    = 1.5d0 * pion*deni
    deiondd = (1.5d0 * dpiondd - eion)*deni
    deiondt = 1.5d0 * dpiondt*deni
    deionda = 1.5d0 * dpionda*deni
    deiondz = 0.0d0


    ! sackur-tetrode equation for the ion entropy of
    ! a single ideal gas characterized by abar
    x       = abar*abar*sqrt(abar) * deni/avo
    s       = sioncon * temp
    z       = x * s * sqrt(s)
    y       = log(z)

    !        y       = 1.0d0/(abar*kt)
    !        yy      = y * sqrt(y)
    !        z       = xni * sifac * yy
    !        etaion  = log(z)


    sion    = (pion*deni + eion)*tempi + kergavo * ytot1 * y
    dsiondd = (dpiondd*deni - pion*deni*deni + deiondd)*tempi &
         - kergavo * deni * ytot1
    dsiondt = (dpiondt*deni + deiondt)*tempi - &
         (pion*deni + eion) * tempi*tempi &
         + 1.5d0 * kergavo * tempi*ytot1
    x       = avo*kerg/abar
    dsionda = (dpionda*deni + deionda)*tempi &
         + kergavo*ytot1*ytot1* (2.5d0 - y)
    dsiondz = 0.0d0



    ! electron-positron section:


    ! assume complete ionization
    xnem    = xni * zbar


    ! enter the table with ye*den
    din = ye*den


    ! bomb proof the input
    temp = min( temp, t(jmax))
    temp = max( temp, t(1))
    din  = min( din , d(imax))
    din  = max( din , d(1))

    ! if (temp .gt. t(jmax)) then
    !    write(6,01) 'temp=',temp,' t(jmax)=',t(jmax)
    !    write(6,*) 'temp too hot, off grid'
    !    write(6,*) 'setting eosfail to true and returning'
    !    eosfail = .true.
    !    return
    ! end if
    ! if (temp .lt. t(1)) then
    !    write(6,01) 'temp=',temp,' t(1)=',t(1)
    !    write(6,*) 'temp too cold, off grid'
    !    write(6,*) 'setting eosfail to true and returning'
    !    eosfail = .true.
    !    return
    ! end if
    ! if (din  .gt. d(imax)) then
    !    write(6,01) 'den*ye=',din,' d(imax)=',d(imax)
    !    write(6,*) 'ye*den too big, off grid'
    !    write(6,*) 'setting eosfail to true and returning'
    !    eosfail = .true.
    !    return
    ! end if
    ! if (din  .lt. d(1)) then
    !    write(6,01) 'ye*den=',din,' d(1)=',d(1)
    !    write(6,*) 'ye*den too small, off grid'
    !    write(6,*) 'setting eosfail to true and returning'
    !    eosfail = .true.
    !    return
    ! end if

    ! hash locate this temperature and density
    jat = int((log10(temp) - tlo)*tstpi) + 1
    jat = max(1,min(jat,jmax-1))
    iat = int((log10(din) - dlo)*dstpi) + 1
    iat = max(1,min(iat,imax-1))


    ! access the table locations only once
    fi(1)  = f(iat,jat)
    fi(2)  = f(iat+1,jat)
    fi(3)  = f(iat,jat+1)
    fi(4)  = f(iat+1,jat+1)
    fi(5)  = ft(iat,jat)
    fi(6)  = ft(iat+1,jat)
    fi(7)  = ft(iat,jat+1)
    fi(8)  = ft(iat+1,jat+1)
    fi(9)  = ftt(iat,jat)
    fi(10) = ftt(iat+1,jat)
    fi(11) = ftt(iat,jat+1)
    fi(12) = ftt(iat+1,jat+1)
    fi(13) = fd(iat,jat)
    fi(14) = fd(iat+1,jat)
    fi(15) = fd(iat,jat+1)
    fi(16) = fd(iat+1,jat+1)
    fi(17) = fdd(iat,jat)
    fi(18) = fdd(iat+1,jat)
    fi(19) = fdd(iat,jat+1)
    fi(20) = fdd(iat+1,jat+1)
    fi(21) = fdt(iat,jat)
    fi(22) = fdt(iat+1,jat)
    fi(23) = fdt(iat,jat+1)
    fi(24) = fdt(iat+1,jat+1)
    fi(25) = fddt(iat,jat)
    fi(26) = fddt(iat+1,jat)
    fi(27) = fddt(iat,jat+1)
    fi(28) = fddt(iat+1,jat+1)
    fi(29) = fdtt(iat,jat)
    fi(30) = fdtt(iat+1,jat)
    fi(31) = fdtt(iat,jat+1)
    fi(32) = fdtt(iat+1,jat+1)
    fi(33) = fddtt(iat,jat)
    fi(34) = fddtt(iat+1,jat)
    fi(35) = fddtt(iat,jat+1)
    fi(36) = fddtt(iat+1,jat+1)


    ! various differences
    xt  = max( (temp - t(jat))*dti_sav(jat), 0.0d0)
    xd  = max( (din - d(iat))*ddi_sav(iat), 0.0d0)
    mxt = 1.0d0 - xt
    mxd = 1.0d0 - xd

    ! the six density and six temperature basis functions
    si0t =   psi0(xt)
    si1t =   psi1(xt)*dt_sav(jat)
    si2t =   psi2(xt)*dt2_sav(jat)

    si0mt =  psi0(mxt)
    si1mt = -psi1(mxt)*dt_sav(jat)
    si2mt =  psi2(mxt)*dt2_sav(jat)

    si0d =   psi0(xd)
    si1d =   psi1(xd)*dd_sav(iat)
    si2d =   psi2(xd)*dd2_sav(iat)

    si0md =  psi0(mxd)
    si1md = -psi1(mxd)*dd_sav(iat)
    si2md =  psi2(mxd)*dd2_sav(iat)

    ! derivatives of the weight functions
    dsi0t =   dpsi0(xt)*dti_sav(jat)
    dsi1t =   dpsi1(xt)
    dsi2t =   dpsi2(xt)*dt_sav(jat)

    dsi0mt = -dpsi0(mxt)*dti_sav(jat)
    dsi1mt =  dpsi1(mxt)
    dsi2mt = -dpsi2(mxt)*dt_sav(jat)

    dsi0d =   dpsi0(xd)*ddi_sav(iat)
    dsi1d =   dpsi1(xd)
    dsi2d =   dpsi2(xd)*dd_sav(iat)

    dsi0md = -dpsi0(mxd)*ddi_sav(iat)
    dsi1md =  dpsi1(mxd)
    dsi2md = -dpsi2(mxd)*dd_sav(iat)

    ! second derivatives of the weight functions
    ddsi0t =   ddpsi0(xt)*dt2i_sav(jat)
    ddsi1t =   ddpsi1(xt)*dti_sav(jat)
    ddsi2t =   ddpsi2(xt)

    ddsi0mt =  ddpsi0(mxt)*dt2i_sav(jat)
    ddsi1mt = -ddpsi1(mxt)*dti_sav(jat)
    ddsi2mt =  ddpsi2(mxt)

    !        ddsi0d =   ddpsi0(xd)*dd2i_sav(iat)
    !        ddsi1d =   ddpsi1(xd)*ddi_sav(iat)
    !        ddsi2d =   ddpsi2(xd)

    !        ddsi0md =  ddpsi0(mxd)*dd2i_sav(iat)
    !        ddsi1md = -ddpsi1(mxd)*ddi_sav(iat)
    !        ddsi2md =  ddpsi2(mxd)


    ! the free energy
    free  = h5(iat,jat, &
         si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
         si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

    ! derivative with respect to density
    df_d  = h5(iat,jat, &
         si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
         dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)


    ! derivative with respect to temperature
    df_t = h5(iat,jat, &
         dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
         si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

    ! derivative with respect to density**2
    !        df_dd = h5(iat,jat,
    !     1          si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
    !     2          ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md)

    ! derivative with respect to temperature**2
    df_tt = h5(iat,jat, &
         ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, &
         si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

    ! derivative with respect to temperature and density
    df_dt = h5(iat,jat, &
         dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
         dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)



    ! now get the pressure derivative with density, chemical potential, and
    ! electron positron number densities
    ! get the interpolation weight functions
    si0t   =  xpsi0(xt)
    si1t   =  xpsi1(xt)*dt_sav(jat)

    si0mt  =  xpsi0(mxt)
    si1mt  =  -xpsi1(mxt)*dt_sav(jat)

    si0d   =  xpsi0(xd)
    si1d   =  xpsi1(xd)*dd_sav(iat)

    si0md  =  xpsi0(mxd)
    si1md  =  -xpsi1(mxd)*dd_sav(iat)


    ! derivatives of weight functions
    dsi0t  = xdpsi0(xt)*dti_sav(jat)
    dsi1t  = xdpsi1(xt)

    dsi0mt = -xdpsi0(mxt)*dti_sav(jat)
    dsi1mt = xdpsi1(mxt)

    dsi0d  = xdpsi0(xd)*ddi_sav(iat)
    dsi1d  = xdpsi1(xd)

    dsi0md = -xdpsi0(mxd)*ddi_sav(iat)
    dsi1md = xdpsi1(mxd)


    ! look in the pressure derivative only once
    fi(1)  = dpdf(iat,jat)
    fi(2)  = dpdf(iat+1,jat)
    fi(3)  = dpdf(iat,jat+1)
    fi(4)  = dpdf(iat+1,jat+1)
    fi(5)  = dpdft(iat,jat)
    fi(6)  = dpdft(iat+1,jat)
    fi(7)  = dpdft(iat,jat+1)
    fi(8)  = dpdft(iat+1,jat+1)
    fi(9)  = dpdfd(iat,jat)
    fi(10) = dpdfd(iat+1,jat)
    fi(11) = dpdfd(iat,jat+1)
    fi(12) = dpdfd(iat+1,jat+1)
    fi(13) = dpdfdt(iat,jat)
    fi(14) = dpdfdt(iat+1,jat)
    fi(15) = dpdfdt(iat,jat+1)
    fi(16) = dpdfdt(iat+1,jat+1)

    ! pressure derivative with density
    dpepdd  = h3(iat,jat, &
         si0t,   si1t,   si0mt,   si1mt, &
         si0d,   si1d,   si0md,   si1md)
    dpepdd  = max(ye * dpepdd,1.0d-30)



    ! look in the electron chemical potential table only once
    fi(1)  = ef(iat,jat)
    fi(2)  = ef(iat+1,jat)
    fi(3)  = ef(iat,jat+1)
    fi(4)  = ef(iat+1,jat+1)
    fi(5)  = eft(iat,jat)
    fi(6)  = eft(iat+1,jat)
    fi(7)  = eft(iat,jat+1)
    fi(8)  = eft(iat+1,jat+1)
    fi(9)  = efd(iat,jat)
    fi(10) = efd(iat+1,jat)
    fi(11) = efd(iat,jat+1)
    fi(12) = efd(iat+1,jat+1)
    fi(13) = efdt(iat,jat)
    fi(14) = efdt(iat+1,jat)
    fi(15) = efdt(iat,jat+1)
    fi(16) = efdt(iat+1,jat+1)


    ! electron chemical potential etaele
    etaele  = h3(iat,jat, &
         si0t,   si1t,   si0mt,   si1mt, &
         si0d,   si1d,   si0md,   si1md)


    ! derivative with respect to density
    x       = h3(iat,jat, &
         si0t,   si1t,   si0mt,   si1mt, &
         dsi0d,  dsi1d,  dsi0md,  dsi1md)
    detadd  = ye * x

    ! derivative with respect to temperature
    detadt  = h3(iat,jat, &
         dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
         si0d,   si1d,   si0md,   si1md)

    ! derivative with respect to abar and zbar
    detada = -x * din * ytot1
    detadz =  x * den * ytot1



    ! look in the number density table only once
    fi(1)  = xf(iat,jat)
    fi(2)  = xf(iat+1,jat)
    fi(3)  = xf(iat,jat+1)
    fi(4)  = xf(iat+1,jat+1)
    fi(5)  = xft(iat,jat)
    fi(6)  = xft(iat+1,jat)
    fi(7)  = xft(iat,jat+1)
    fi(8)  = xft(iat+1,jat+1)
    fi(9)  = xfd(iat,jat)
    fi(10) = xfd(iat+1,jat)
    fi(11) = xfd(iat,jat+1)
    fi(12) = xfd(iat+1,jat+1)
    fi(13) = xfdt(iat,jat)
    fi(14) = xfdt(iat+1,jat)
    fi(15) = xfdt(iat,jat+1)
    fi(16) = xfdt(iat+1,jat+1)

    ! electron + positron number densities
    xnefer   = h3(iat,jat, &
         si0t,   si1t,   si0mt,   si1mt, &
         si0d,   si1d,   si0md,   si1md)

    ! derivative with respect to density
    x        = h3(iat,jat, &
         si0t,   si1t,   si0mt,   si1mt, &
         dsi0d,  dsi1d,  dsi0md,  dsi1md)
    x = max(x,1.0d-30)
    dxnedd   = ye * x

    ! derivative with respect to temperature
    dxnedt   = h3(iat,jat, &
         dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
         si0d,   si1d,   si0md,   si1md)

    ! derivative with respect to abar and zbar
    dxneda = -x * din * ytot1
    dxnedz =  x  * den * ytot1


    ! the desired electron-positron thermodynamic quantities

    ! dpepdd at high temperatures and low densities is below the
    ! floating point limit of the subtraction of two large terms.
    ! since dpresdd doesn't enter the maxwell relations at all, use the
    ! bicubic interpolation done above instead of the formally correct expression
    x       = din * din
    pele    = x * df_d
    dpepdt  = x * df_dt
    !        dpepdd  = ye * (x * df_dd + 2.0d0 * din * df_d)
    s       = dpepdd/ye - 2.0d0 * din * df_d
    dpepda  = -ytot1 * (2.0d0 * pele + s * din)
    dpepdz  = den*ytot1*(2.0d0 * din * df_d  +  s)


    x       = ye * ye
    sele    = -df_t * ye
    dsepdt  = -df_tt * ye
    dsepdd  = -df_dt * x
    dsepda  = ytot1 * (ye * df_dt * din - sele)
    dsepdz  = -ytot1 * (ye * df_dt * den  + df_t)


    eele    = ye*free + temp * sele
    deepdt  = temp * dsepdt
    deepdd  = x * df_d + temp * dsepdd
    deepda  = -ye * ytot1 * (free +  df_d * din) + temp * dsepda
    deepdz  = ytot1* (free + ye * df_d * den) + temp * dsepdz




    ! coulomb section:

    ! uniform background corrections only
    ! from yakovlev & shalybkov 1989
    ! lami is the average ion seperation
    ! plasg is the plasma coupling parameter

    z        = forth * pi
    s        = z * xni
    dsdd     = z * dxnidd
    dsda     = z * dxnida

    lami     = 1.0d0/s**third
    inv_lami = 1.0d0/lami
    z        = -third * lami
    lamidd   = z * dsdd/s
    lamida   = z * dsda/s

    plasg    = zbar*zbar*esqu*ktinv*inv_lami
    z        = -plasg * inv_lami
    plasgdd  = z * lamidd
    plasgda  = z * lamida
    plasgdt  = -plasg*ktinv * kerg
    plasgdz  = 2.0d0 * plasg/zbar


    ! ! yakovlev & shalybkov 1989 equations 82, 85, 86, 87
    ! if (plasg .ge. 1.0) then
    !    x        = plasg**(0.25d0)
    !    y        = avo * ytot1 * kerg
    !    ecoul    = y * temp * (a1*plasg + b1*x + c1/x + d1)
    !    pcoul    = third * den * ecoul
    !    scoul    = -y * (3.0d0*b1*x - 5.0d0*c1/x &
    !         + d1 * (log(plasg) - 1.0d0) - e1)

    !    y        = avo*ytot1*kt*(a1 + 0.25d0/plasg*(b1*x - c1/x))
    !    decouldd = y * plasgdd
    !    decouldt = y * plasgdt + ecoul/temp
    !    decoulda = y * plasgda - ecoul/abar
    !    decouldz = y * plasgdz

    !    y        = third * den
    !    dpcouldd = third * ecoul + y*decouldd
    !    dpcouldt = y * decouldt
    !    dpcoulda = y * decoulda
    !    dpcouldz = y * decouldz


    !    y        = -avo*kerg/(abar*plasg)*(0.75d0*b1*x+1.25d0*c1/x+d1)
    !    dscouldd = y * plasgdd
    !    dscouldt = y * plasgdt
    !    dscoulda = y * plasgda - scoul/abar
    !    dscouldz = y * plasgdz


    !    ! yakovlev & shalybkov 1989 equations 102, 103, 104
    ! else if (plasg .lt. 1.0) then
    !    x        = plasg*sqrt(plasg)
    !    y        = plasg**b2
    !    z        = c2 * x - third * a2 * y
    !    pcoul    = -pion * z
    !    ecoul    = 3.0d0 * pcoul/den
    !    scoul    = -avo/abar*kerg*(c2*x -a2*(b2-1.0d0)/b2*y)

    !    s        = 1.5d0*c2*x/plasg - third*a2*b2*y/plasg
    !    dpcouldd = -dpiondd*z - pion*s*plasgdd
    !    dpcouldt = -dpiondt*z - pion*s*plasgdt
    !    dpcoulda = -dpionda*z - pion*s*plasgda
    !    dpcouldz = -dpiondz*z - pion*s*plasgdz

    !    s        = 3.0d0/den
    !    decouldd = s * dpcouldd - ecoul/den
    !    decouldt = s * dpcouldt
    !    decoulda = s * dpcoulda
    !    decouldz = s * dpcouldz

    !    s        = -avo*kerg/(abar*plasg)*(1.5d0*c2*x-a2*(b2-1.0d0)*y)
    !    dscouldd = s * plasgdd
    !    dscouldt = s * plasgdt
    !    dscoulda = s * plasgda - scoul/abar
    !    dscouldz = s * plasgdz
    ! end if


    ! ! bomb proof
    ! x   = prad + pion + pele + pcoul
    ! y   = erad + eion + eele + ecoul
    ! z   = srad + sion + sele + scoul

    ! !        write(6,*) x,y,z
    ! !        if (x .le. 0.0 .or. y .le. 0.0 .or. z .le. 0.0) then
    ! if (x .le. 0.0 .or. y .le. 0.0) then
    !    !        if (x .le. 0.0) then

    !    !         write(6,*)
    !    !         write(6,*) 'coulomb corrections are causing a negative pressure'
    !    !         write(6,*) 'setting all coulomb corrections to zero'
    !    !         write(6,*)

    !    pcoul    = 0.0d0
    !    dpcouldd = 0.0d0
    !    dpcouldt = 0.0d0
    !    dpcoulda = 0.0d0
    !    dpcouldz = 0.0d0
    !    ecoul    = 0.0d0
    !    decouldd = 0.0d0
    !    decouldt = 0.0d0
    !    decoulda = 0.0d0
    !    decouldz = 0.0d0
    !    scoul    = 0.0d0
    !    dscouldd = 0.0d0
    !    dscouldt = 0.0d0
    !    dscoulda = 0.0d0
    !    dscouldz = 0.0d0
    ! end if

    pcoul    = 0.0d0
    dpcouldd = 0.0d0
    dpcouldt = 0.0d0
    dpcoulda = 0.0d0
    dpcouldz = 0.0d0
    ecoul    = 0.0d0
    decouldd = 0.0d0
    decouldt = 0.0d0
    decoulda = 0.0d0
    decouldz = 0.0d0
    scoul    = 0.0d0
    dscouldd = 0.0d0
    dscouldt = 0.0d0
    dscoulda = 0.0d0
    dscouldz = 0.0d0

    ! sum all the gas components
    pgas    = pion + pele + pcoul
    egas    = eion + eele + ecoul
    sgas    = sion + sele + scoul

    dpgasdd = dpiondd + dpepdd + dpcouldd
    dpgasdt = dpiondt + dpepdt + dpcouldt
    dpgasda = dpionda + dpepda + dpcoulda
    dpgasdz = dpiondz + dpepdz + dpcouldz

    degasdd = deiondd + deepdd + decouldd
    degasdt = deiondt + deepdt + decouldt
    degasda = deionda + deepda + decoulda
    degasdz = deiondz + deepdz + decouldz

    dsgasdd = dsiondd + dsepdd + dscouldd
    dsgasdt = dsiondt + dsepdt + dscouldt
    dsgasda = dsionda + dsepda + dscoulda
    dsgasdz = dsiondz + dsepdz + dscouldz




    ! add in radiation to get the total
    pres    = prad + pgas
    ener    = erad + egas
    entr    = srad + sgas

    dpresdd = dpraddd + dpgasdd
    dpresdt = dpraddt + dpgasdt
    dpresda = dpradda + dpgasda
    dpresdz = dpraddz + dpgasdz

    denerdd = deraddd + degasdd
    denerdt = deraddt + degasdt
    denerda = deradda + degasda
    denerdz = deraddz + degasdz

    dentrdd = dsraddd + dsgasdd
    dentrdt = dsraddt + dsgasdt
    dentrda = dsradda + dsgasda
    dentrdz = dsraddz + dsgasdz


    ! for the gas
    ! the temperature and density exponents (c&g 9.81 9.82)
    ! the specific heat at constant volume (c&g 9.92)
    ! the third adiabatic exponent (c&g 9.93)
    ! the first adiabatic exponent (c&g 9.97)
    ! the second adiabatic exponent (c&g 9.105)
    ! the specific heat at constant pressure (c&g 9.98)
    ! and relativistic formula for the sound speed (c&g 14.29)

    zz        = pgas*deni
    zzi       = den/pgas
    chit_gas  = temp/pgas * dpgasdt
    chid_gas  = dpgasdd*zzi
    cv_gas    = degasdt
    x         = zz * chit_gas/(temp * cv_gas)
    gam3_gas  = x + 1.0d0
    gam1_gas  = chit_gas*x + chid_gas
    nabad_gas = x/gam1_gas
    gam2_gas  = 1.0d0/(1.0d0 - nabad_gas)
    cp_gas    = cv_gas * gam1_gas/chid_gas
    z         = 1.0d0 + (egas + light2)*zzi
    sound_gas = clight * sqrt(gam1_gas/z)



    ! for the totals
    zz    = pres*deni
    zzi   = den/pres
    chit  = temp/pres * dpresdt
    chid  = dpresdd*zzi
    cv    = denerdt
    x     = zz * chit/(temp * cv)
    gam3  = x + 1.0d0
    gam1  = chit*x + chid
    nabad = x/gam1
    gam2  = 1.0d0/(1.0d0 - nabad)
    cp    = cv * gam1/chid
    ! modified from the original helmeos
    eps   = ener + mexc_in*mev_to_erg/m_u_g
    z     = 1.0d0 + (eps + light2)*zzi
    sound = clight * sqrt(gam1/z)
    
    if(present(chi_out))then
       chi_out = dpresdd - dpresdt*denerdd/denerdt
       ! chi_out = pres*deni * (gam1 - x)
    endif
    
    if(present(kappa_out))then
       kappa_out = dpresdt / (den*denerdt)
       ! kappa_out = x
    endif

    ! maxwell relations; each is zero if the consistency is perfect
    x   = den * den

    dse = temp*dentrdt/denerdt - 1.0d0

    dpe = (denerdd*x + temp*dpresdt)/pres - 1.0d0

    dsp = -dentrdd*x/dpresdt - 1.0d0
    
    
    eps_out = ener + mexc_in*mev_to_erg/m_u_g
    pres_out = pres
    cs2_out = sound**2
    entr_out = entr * amu/kerg

    
    ! ! store this row
    ! ptot_row(j)   = pres
    ! dpt_row(j)    = dpresdt
    ! dpd_row(j)    = dpresdd
    ! dpa_row(j)    = dpresda
    ! dpz_row(j)    = dpresdz
    
    ! etot_row(j)   = ener
    ! det_row(j)    = denerdt
    ! ded_row(j)    = denerdd
    ! dea_row(j)    = denerda
    ! dez_row(j)    = denerdz
    
    ! stot_row(j)   = entr
    ! dst_row(j)    = dentrdt
    ! dsd_row(j)    = dentrdd
    ! dsa_row(j)    = dentrda
    ! dsz_row(j)    = dentrdz
    
    
    ! pgas_row(j)   = pgas
    ! dpgast_row(j) = dpgasdt
    ! dpgasd_row(j) = dpgasdd
    ! dpgasa_row(j) = dpgasda
    ! dpgasz_row(j) = dpgasdz
    
    ! egas_row(j)   = egas
    ! degast_row(j) = degasdt
    ! degasd_row(j) = degasdd
    ! degasa_row(j) = degasda
    ! degasz_row(j) = degasdz
    
    ! sgas_row(j)   = sgas
    ! dsgast_row(j) = dsgasdt
    ! dsgasd_row(j) = dsgasdd
    ! dsgasa_row(j) = dsgasda
    ! dsgasz_row(j) = dsgasdz
    
    
    ! prad_row(j)   = prad
    ! dpradt_row(j) = dpraddt
    ! dpradd_row(j) = dpraddd
    ! dprada_row(j) = dpradda
    ! dpradz_row(j) = dpraddz
    
    ! erad_row(j)   = erad
    ! deradt_row(j) = deraddt
    ! deradd_row(j) = deraddd
    ! derada_row(j) = deradda
    ! deradz_row(j) = deraddz
    
    ! srad_row(j)   = srad
    ! dsradt_row(j) = dsraddt
    ! dsradd_row(j) = dsraddd
    ! dsrada_row(j) = dsradda
    ! dsradz_row(j) = dsraddz
    
    
    ! pion_row(j)   = pion
    ! dpiont_row(j) = dpiondt
    ! dpiond_row(j) = dpiondd
    ! dpiona_row(j) = dpionda
    ! dpionz_row(j) = dpiondz
    
    ! eion_row(j)   = eion
    ! deiont_row(j) = deiondt
    ! deiond_row(j) = deiondd
    ! deiona_row(j) = deionda
    ! deionz_row(j) = deiondz
    
    ! sion_row(j)   = sion
    ! dsiont_row(j) = dsiondt
    ! dsiond_row(j) = dsiondd
    ! dsiona_row(j) = dsionda
    ! dsionz_row(j) = dsiondz
    
    ! xni_row(j)    = xni
    
    ! pele_row(j)   = pele
    ! ppos_row(j)   = 0.0d0
    ! dpept_row(j)  = dpepdt
    ! dpepd_row(j)  = dpepdd
    ! dpepa_row(j)  = dpepda
    ! dpepz_row(j)  = dpepdz
    
    ! eele_row(j)   = eele
    ! epos_row(j)   = 0.0d0
    ! deept_row(j)  = deepdt
    ! deepd_row(j)  = deepdd
    ! deepa_row(j)  = deepda
    ! deepz_row(j)  = deepdz
    
    ! sele_row(j)   = sele
    ! spos_row(j)   = 0.0d0
    ! dsept_row(j)  = dsepdt
    ! dsepd_row(j)  = dsepdd
    ! dsepa_row(j)  = dsepda
    ! dsepz_row(j)  = dsepdz
    
    ! xnem_row(j)   = xnem
    ! xne_row(j)    = xnefer
    ! dxnet_row(j)  = dxnedt
    ! dxned_row(j)  = dxnedd
    ! dxnea_row(j)  = dxneda
    ! dxnez_row(j)  = dxnedz
    ! xnp_row(j)    = 0.0d0
    ! zeff_row(j)   = zbar
    
    ! etaele_row(j) = etaele
    ! detat_row(j)  = detadt
    ! detad_row(j)  = detadd
    ! detaa_row(j)  = detada
    ! detaz_row(j)  = detadz
    ! etapos_row(j) = 0.0d0
    
    ! pcou_row(j)   = pcoul
    ! dpcout_row(j) = dpcouldt
    ! dpcoud_row(j) = dpcouldd
    ! dpcoua_row(j) = dpcoulda
    ! dpcouz_row(j) = dpcouldz
    
    ! ecou_row(j)   = ecoul
    ! decout_row(j) = decouldt
    ! decoud_row(j) = decouldd
    ! decoua_row(j) = decoulda
    ! decouz_row(j) = decouldz
    
    ! scou_row(j)   = scoul
    ! dscout_row(j) = dscouldt
    ! dscoud_row(j) = dscouldd
    ! dscoua_row(j) = dscoulda
    ! dscouz_row(j) = dscouldz
    
    ! plasg_row(j)  = plasg
    
    ! dse_row(j)    = dse
    ! dpe_row(j)    = dpe
    ! dsp_row(j)    = dsp
    
    ! cv_gas_row(j)    = cv_gas
    ! cp_gas_row(j)    = cp_gas
    ! gam1_gas_row(j)  = gam1_gas
    ! gam2_gas_row(j)  = gam2_gas
    ! gam3_gas_row(j)  = gam3_gas
    ! nabad_gas_row(j) = nabad_gas
    ! cs_gas_row(j)    = sound_gas
    
    ! cv_row(j)     = cv
    ! cp_row(j)     = cp
    ! gam1_row(j)   = gam1
    ! gam2_row(j)   = gam2
    ! gam3_row(j)   = gam3
    ! nabad_row(j)  = nabad
    ! cs_row(j)     = sound
    
  end subroutine eos_all


  subroutine eos_get_hmin(mexc, hmin)
    real(8),intent(in) :: mexc
    real(8),intent(out) :: hmin

    hmin = clight**2 * (1d0 + mexc/m_u_mev)
    
  end subroutine eos_get_hmin


  subroutine eos_get_eps_min(rho, ye, mexc, eps_min)
    real(8),intent(in) :: rho, ye, mexc
    real(8),intent(out) :: eps_min

    real(8) :: ne,x,chi,eps_e

    ne = rho*ye/m_u_g
    x = lambda_e*(3d0*pi**2*ne)**(1d0/3d0)

    chi = ( x*sqrt(1d0+x**2)*(2d0/3d0*x**2-1d0) + log(x+sqrt(1d0+x**2)) ) / (8d0*pi**2)
    
    eps_e = m_e_erg/lambda_e3 * chi / rho

    eps_min = eps_e + mexc/m_u_mev*clight**2
    
  end subroutine eos_get_eps_min


  subroutine eos_get_range(temp_min_out,temp_max_out,rhoye_min_out,rhoye_max_out)
    real(8),intent(out) :: temp_min_out,temp_max_out,rhoye_min_out,rhoye_max_out

    temp_min_out = temp_eos_min
    temp_max_out = temp_eos_max
    rhoye_min_out = d(1)
    rhoye_max_out = d(imax)
    
  end subroutine eos_get_range


  subroutine read_helm_table
    ! include 'implno.dek'
    ! include 'helm_table_storage.dek'

    ! this routine reads the helmholtz eos file, and
    ! must be called once before the helmeos routine is invoked.

    ! declare local variables
    integer          i,j
    double precision tsav,dsav,dth,dt2,dti,dt2i,dt3i, &
         dd,dd2,ddi,dd2i,dd3i


    ! open the file (use softlinks to input the desired table)

    open(unit=19,file=fn_eos,status='old')

    ! for standard table limits
    tlo   = 3.0d0
    thi   = 13.0d0
    tstp  = (thi - tlo)/dble(jmax-1)
    tstpi = 1.0d0/tstp
    dlo   = -12.0d0
    dhi   = 15.0d0
    dstp  = (dhi - dlo)/dble(imax-1)
    dstpi = 1.0d0/dstp

    ! read the helmholtz free energy and its derivatives
    do j=1,jmax
       tsav = tlo + (j-1)*tstp
       t(j) = 10.0d0**(tsav)
       do i=1,imax
          dsav = dlo + (i-1)*dstp
          d(i) = 10.0d0**(dsav)
          read(19,*) f(i,j),fd(i,j),ft(i,j),fdd(i,j),ftt(i,j),fdt(i,j), &
               fddt(i,j),fdtt(i,j),fddtt(i,j)
       enddo
    enddo
    !       write(6,*) 'read main table'


    ! read the pressure derivative with density table
    do j=1,jmax
       do i=1,imax
          read(19,*) dpdf(i,j),dpdfd(i,j),dpdft(i,j),dpdfdt(i,j)
       enddo
    enddo
    !       write(6,*) 'read dpdd table'

    ! read the electron chemical potential table
    do j=1,jmax
       do i=1,imax
          read(19,*) ef(i,j),efd(i,j),eft(i,j),efdt(i,j)
       enddo
    enddo
    !       write(6,*) 'read eta table'

    ! read the number density table
    do j=1,jmax
       do i=1,imax
          read(19,*) xf(i,j),xfd(i,j),xft(i,j),xfdt(i,j)
       enddo
    enddo
    !       write(6,*) 'read xne table'

    ! close the file
    close(unit=19)


    ! construct the temperature and density deltas and their inverses
    do j=1,jmax-1
       dth          = t(j+1) - t(j)
       dt2         = dth * dth
       dti         = 1.0d0/dth
       dt2i        = 1.0d0/dt2
       dt3i        = dt2i*dti
       dt_sav(j)   = dth
       dt2_sav(j)  = dt2
       dti_sav(j)  = dti
       dt2i_sav(j) = dt2i
       dt3i_sav(j) = dt3i
    end do
    do i=1,imax-1
       dd          = d(i+1) - d(i)
       dd2         = dd * dd
       ddi         = 1.0d0/dd
       dd2i        = 1.0d0/dd2
       dd3i        = dd2i*ddi
       dd_sav(i)   = dd
       dd2_sav(i)  = dd2
       ddi_sav(i)  = ddi
       dd2i_sav(i) = dd2i
       dd3i_sav(i) = dd3i
    enddo



    !      write(6,*)
    !      write(6,*) 'finished reading eos table'
    !      write(6,04) 'imax=',imax,' jmax=',jmax
    !04    format(1x,4(a,i4))
    !      write(6,03) 'temp(1)   =',t(1),' temp(jmax)   =',t(jmax)
    !      write(6,03) 'ye*den(1) =',d(1),' ye*den(imax) =',d(imax)
    !03    format(1x,4(a,1pe11.3))
    !      write(6,*)

    return
  end subroutine read_helm_table



  subroutine helmeos
    !include 'implno.dek'
    !include 'const.dek'
    !include 'vector_eos.dek'
    !include 'helm_table_storage.dek'


    ! given a temperature temp [K], density den [g/cm**3], and a composition
    ! characterized by abar and zbar, this routine returns most of the other
    ! thermodynamic quantities. of prime interest is the pressure [erg/cm**3],
    ! specific thermal energy [erg/gr], the entropy [erg/g/K], along with
    ! their derivatives with respect to temperature, density, abar, and zbar.
    ! other quantites such the normalized chemical potential eta (plus its
    ! derivatives), number density of electrons and positron pair (along
    ! with their derivatives), adiabatic indices, specific heats, and
    ! relativistically correct sound speed are also returned.
    !
    ! this routine assumes planckian photons, an ideal gas of ions,
    ! and an electron-positron gas with an arbitrary degree of relativity
    ! and degeneracy. interpolation in a table of the helmholtz free energy
    ! is used to return the electron-positron thermodynamic quantities.
    ! all other derivatives are analytic.
    !
    ! references: cox & giuli chapter 24 ; timmes & swesty apj 1999


    ! declare
    integer          i,j
    double precision temp,den,abar,zbar,ytot1,ye, &
         x,y,zz,zzi,deni,tempi,xni,dxnidd,dxnida, &
         dpepdt,dpepdd,deepdt,deepdd,dsepdd,dsepdt, &
         dpraddd,dpraddt,deraddd,deraddt,dpiondd,dpiondt, &
         deiondd,deiondt,dsraddd,dsraddt,dsiondd,dsiondt, &
         dse,dpe,dsp,kt,ktinv,prad,erad,srad,pion,eion, &
         sion,xnem,pele,eele,sele,pres,ener,entr,dpresdd, &
         dpresdt,denerdd,denerdt,dentrdd,dentrdt,cv,cp, &
         gam1,gam2,gam3,chit,chid,nabad,sound,etaele, &
         detadt,detadd,xnefer,dxnedt,dxnedd,s

    double precision pgas,dpgasdd,dpgasdt,dpgasda,dpgasdz, &
         egas,degasdd,degasdt,degasda,degasdz, &
         sgas,dsgasdd,dsgasdt,dsgasda,dsgasdz, &
         cv_gas,cp_gas,gam1_gas,gam2_gas,gam3_gas, &
         chit_gas,chid_gas,nabad_gas,sound_gas


    double precision sioncon,forth,forpi,kergavo,ikavo,asoli3,light2
    parameter        (sioncon = (2.0d0 * pi * amu * kerg)/(h*h), &
         forth   = 4.0d0/3.0d0, &
         forpi   = 4.0d0 * pi, &
         kergavo = kerg * avo, &
         ikavo   = 1.0d0/kergavo, &
         asoli3  = asol/3.0d0, &
         light2  = clight * clight)

    ! for the abar derivatives
    double precision dpradda,deradda,dsradda, &
         dpionda,deionda,dsionda, &
         dpepda,deepda,dsepda, &
         dpresda,denerda,dentrda, &
         detada,dxneda

    ! for the zbar derivatives
    double precision dpraddz,deraddz,dsraddz, &
         dpiondz,deiondz,dsiondz, &
         dpepdz,deepdz,dsepdz, &
         dpresdz,denerdz,dentrdz, &
         detadz,dxnedz

    ! for the interpolations
    integer          iat,jat
    double precision free,df_d,df_t,df_dd,df_tt,df_dt
    double precision xt,xd,mxt,mxd, &
         si0t,si1t,si2t,si0mt,si1mt,si2mt, &
         si0d,si1d,si2d,si0md,si1md,si2md, &
         dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt, &
         dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md, &
         ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt, &
         ddsi0d,ddsi1d,ddsi2d,ddsi0md,ddsi1md,ddsi2md, &
         z,psi0,dpsi0,ddpsi0,psi1,dpsi1,ddpsi1,psi2, &
         dpsi2,ddpsi2,din,h5,fi(36), &
         xpsi0,xdpsi0,xpsi1,xdpsi1,h3, &
         w0t,w1t,w2t,w0mt,w1mt,w2mt, &
         w0d,w1d,w2d,w0md,w1md,w2md


    ! for the uniform background coulomb correction
    double precision dsdd,dsda,lami,inv_lami,lamida,lamidd, &
         plasg,plasgdd,plasgdt,plasgda,plasgdz, &
         ecoul,decouldd,decouldt,decoulda,decouldz, &
         pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz, &
         scoul,dscouldd,dscouldt,dscoulda,dscouldz, &
         a1,b1,c1,d1,e1,a2,b2,c2,third,esqu
    parameter        (a1    = -0.898004d0, &
         b1    =  0.96786d0, &
         c1    =  0.220703d0, &
         d1    = -0.86097d0, &
         e1    =  2.5269d0, &
         a2    =  0.29561d0, &
         b2    =  1.9885d0, &
         c2    =  0.288675d0, &
         third =  1.0d0/3.0d0, &
         esqu  =  qe * qe)


    ! quintic hermite polynomial statement functions
    ! psi0 and its derivatives
    psi0(z)   = z**3 * ( z * (-6.0d0*z + 15.0d0) -10.0d0) + 1.0d0
    dpsi0(z)  = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
    ddpsi0(z) = z* ( z*( -120.0d0*z + 180.0d0) -60.0d0)


    ! psi1 and its derivatives
    psi1(z)   = z* ( z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
    dpsi1(z)  = z*z * ( z * (-15.0d0*z + 32.0d0) - 18.0d0) +1.0d0
    ddpsi1(z) = z * (z * (-60.0d0*z + 96.0d0) -36.0d0)


    ! psi2  and its derivatives
    psi2(z)   = 0.5d0*z*z*( z* ( z * (-z + 3.0d0) - 3.0d0) + 1.0d0)
    dpsi2(z)  = 0.5d0*z*( z*(z*(-5.0d0*z + 12.0d0) - 9.0d0) + 2.0d0)
    ddpsi2(z) = 0.5d0*(z*( z * (-20.0d0*z + 36.0d0) - 18.0d0) + 2.0d0)


    ! biquintic hermite polynomial statement function
    h5(i,j,w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md)= &
         fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t &
         + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt &
         + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t &
         + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt &
         + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t &
         + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt &
         + fi(13) *w1d*w0t   + fi(14) *w1md*w0t &
         + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt &
         + fi(17) *w2d*w0t   + fi(18) *w2md*w0t &
         + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt &
         + fi(21) *w1d*w1t   + fi(22) *w1md*w1t &
         + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt &
         + fi(25) *w2d*w1t   + fi(26) *w2md*w1t &
         + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt &
         + fi(29) *w1d*w2t   + fi(30) *w1md*w2t &
         + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt &
         + fi(33) *w2d*w2t   + fi(34) *w2md*w2t &
         + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt



    ! cubic hermite polynomial statement functions
    ! psi0 & derivatives
    xpsi0(z)  = z * z * (2.0d0*z - 3.0d0) + 1.0
    xdpsi0(z) = z * (6.0d0*z - 6.0d0)


    ! psi1 & derivatives
    xpsi1(z)  = z * ( z * (z - 2.0d0) + 1.0d0)
    xdpsi1(z) = z * (3.0d0*z - 4.0d0) + 1.0d0


    ! bicubic hermite polynomial statement function
    h3(i,j,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) = &
         fi(1)  *w0d*w0t   +  fi(2)  *w0md*w0t &
         + fi(3)  *w0d*w0mt  +  fi(4)  *w0md*w0mt &
         + fi(5)  *w0d*w1t   +  fi(6)  *w0md*w1t &
         + fi(7)  *w0d*w1mt  +  fi(8)  *w0md*w1mt &
         + fi(9)  *w1d*w0t   +  fi(10) *w1md*w0t &
         + fi(11) *w1d*w0mt  +  fi(12) *w1md*w0mt &
         + fi(13) *w1d*w1t   +  fi(14) *w1md*w1t &
         + fi(15) *w1d*w1mt  +  fi(16) *w1md*w1mt



    ! popular format statements
01  format(1x,5(a,1pe11.3))
02  format(1x,a,1p4e16.8)
03  format(1x,4(a,1pe11.3))
04  format(1x,4(a,i4))



    ! start of pipeline loop, normal execution starts here
    !eosfail = .false.
    do j=jlo_eos,jhi_eos

       !       if (temp_row(j) .le. 0.0) stop 'temp less than 0 in helmeos'
       !       if (den_row(j)  .le. 0.0) stop 'den less than 0 in helmeos'

       temp  = temp_row(j)
       den   = den_row(j)
       abar  = abar_row(j)
       zbar  = zbar_row(j)
       ytot1 = 1.0d0/abar
       ye    = max(1.0d-16,ytot1 * zbar)



       ! initialize
       deni    = 1.0d0/den
       tempi   = 1.0d0/temp
       kt      = kerg * temp
       ktinv   = 1.0d0/kt


       ! radiation section:
       prad    = asoli3 * temp * temp * temp * temp
       dpraddd = 0.0d0
       dpraddt = 4.0d0 * prad*tempi
       dpradda = 0.0d0
       dpraddz = 0.0d0

       erad    = 3.0d0 * prad*deni
       deraddd = -erad*deni
       deraddt = 3.0d0 * dpraddt*deni
       deradda = 0.0d0
       deraddz = 0.0d0

       srad    = (prad*deni + erad)*tempi
       dsraddd = (dpraddd*deni - prad*deni*deni + deraddd)*tempi
       dsraddt = (dpraddt*deni + deraddt - srad)*tempi
       dsradda = 0.0d0
       dsraddz = 0.0d0


       ! ion section:
       xni     = avo * ytot1 * den
       dxnidd  = avo * ytot1
       dxnida  = -xni * ytot1

       pion    = xni * kt
       dpiondd = dxnidd * kt
       dpiondt = xni * kerg
       dpionda = dxnida * kt
       dpiondz = 0.0d0

       eion    = 1.5d0 * pion*deni
       deiondd = (1.5d0 * dpiondd - eion)*deni
       deiondt = 1.5d0 * dpiondt*deni
       deionda = 1.5d0 * dpionda*deni
       deiondz = 0.0d0


       ! sackur-tetrode equation for the ion entropy of
       ! a single ideal gas characterized by abar
       x       = abar*abar*sqrt(abar) * deni/avo
       s       = sioncon * temp
       z       = x * s * sqrt(s)
       y       = log(z)

       !        y       = 1.0d0/(abar*kt)
       !        yy      = y * sqrt(y)
       !        z       = xni * sifac * yy
       !        etaion  = log(z)


       sion    = (pion*deni + eion)*tempi + kergavo * ytot1 * y
       dsiondd = (dpiondd*deni - pion*deni*deni + deiondd)*tempi &
            - kergavo * deni * ytot1
       dsiondt = (dpiondt*deni + deiondt)*tempi - &
            (pion*deni + eion) * tempi*tempi &
            + 1.5d0 * kergavo * tempi*ytot1
       x       = avo*kerg/abar
       dsionda = (dpionda*deni + deionda)*tempi &
            + kergavo*ytot1*ytot1* (2.5d0 - y)
       dsiondz = 0.0d0



       ! electron-positron section:


       ! assume complete ionization
       xnem    = xni * zbar


       ! enter the table with ye*den
       din = ye*den


       ! ! bomb proof the input
       ! if (temp .gt. t(jmax)) then
       !    write(6,01) 'temp=',temp,' t(jmax)=',t(jmax)
       !    write(6,*) 'temp too hot, off grid'
       !    write(6,*) 'setting eosfail to true and returning'
       !    eosfail = .true.
       !    return
       ! end if
       ! if (temp .lt. t(1)) then
       !    write(6,01) 'temp=',temp,' t(1)=',t(1)
       !    write(6,*) 'temp too cold, off grid'
       !    write(6,*) 'setting eosfail to true and returning'
       !    eosfail = .true.
       !    return
       ! end if
       ! if (din  .gt. d(imax)) then
       !    write(6,01) 'den*ye=',din,' d(imax)=',d(imax)
       !    write(6,*) 'ye*den too big, off grid'
       !    write(6,*) 'setting eosfail to true and returning'
       !    eosfail = .true.
       !    return
       ! end if
       ! if (din  .lt. d(1)) then
       !    write(6,01) 'ye*den=',din,' d(1)=',d(1)
       !    write(6,*) 'ye*den too small, off grid'
       !    write(6,*) 'setting eosfail to true and returning'
       !    eosfail = .true.
       !    return
       ! end if

       ! hash locate this temperature and density
       jat = int((log10(temp) - tlo)*tstpi) + 1
       jat = max(1,min(jat,jmax-1))
       iat = int((log10(din) - dlo)*dstpi) + 1
       iat = max(1,min(iat,imax-1))


       ! access the table locations only once
       fi(1)  = f(iat,jat)
       fi(2)  = f(iat+1,jat)
       fi(3)  = f(iat,jat+1)
       fi(4)  = f(iat+1,jat+1)
       fi(5)  = ft(iat,jat)
       fi(6)  = ft(iat+1,jat)
       fi(7)  = ft(iat,jat+1)
       fi(8)  = ft(iat+1,jat+1)
       fi(9)  = ftt(iat,jat)
       fi(10) = ftt(iat+1,jat)
       fi(11) = ftt(iat,jat+1)
       fi(12) = ftt(iat+1,jat+1)
       fi(13) = fd(iat,jat)
       fi(14) = fd(iat+1,jat)
       fi(15) = fd(iat,jat+1)
       fi(16) = fd(iat+1,jat+1)
       fi(17) = fdd(iat,jat)
       fi(18) = fdd(iat+1,jat)
       fi(19) = fdd(iat,jat+1)
       fi(20) = fdd(iat+1,jat+1)
       fi(21) = fdt(iat,jat)
       fi(22) = fdt(iat+1,jat)
       fi(23) = fdt(iat,jat+1)
       fi(24) = fdt(iat+1,jat+1)
       fi(25) = fddt(iat,jat)
       fi(26) = fddt(iat+1,jat)
       fi(27) = fddt(iat,jat+1)
       fi(28) = fddt(iat+1,jat+1)
       fi(29) = fdtt(iat,jat)
       fi(30) = fdtt(iat+1,jat)
       fi(31) = fdtt(iat,jat+1)
       fi(32) = fdtt(iat+1,jat+1)
       fi(33) = fddtt(iat,jat)
       fi(34) = fddtt(iat+1,jat)
       fi(35) = fddtt(iat,jat+1)
       fi(36) = fddtt(iat+1,jat+1)


       ! various differences
       xt  = max( (temp - t(jat))*dti_sav(jat), 0.0d0)
       xd  = max( (din - d(iat))*ddi_sav(iat), 0.0d0)
       mxt = 1.0d0 - xt
       mxd = 1.0d0 - xd

       ! the six density and six temperature basis functions
       si0t =   psi0(xt)
       si1t =   psi1(xt)*dt_sav(jat)
       si2t =   psi2(xt)*dt2_sav(jat)

       si0mt =  psi0(mxt)
       si1mt = -psi1(mxt)*dt_sav(jat)
       si2mt =  psi2(mxt)*dt2_sav(jat)

       si0d =   psi0(xd)
       si1d =   psi1(xd)*dd_sav(iat)
       si2d =   psi2(xd)*dd2_sav(iat)

       si0md =  psi0(mxd)
       si1md = -psi1(mxd)*dd_sav(iat)
       si2md =  psi2(mxd)*dd2_sav(iat)

       ! derivatives of the weight functions
       dsi0t =   dpsi0(xt)*dti_sav(jat)
       dsi1t =   dpsi1(xt)
       dsi2t =   dpsi2(xt)*dt_sav(jat)

       dsi0mt = -dpsi0(mxt)*dti_sav(jat)
       dsi1mt =  dpsi1(mxt)
       dsi2mt = -dpsi2(mxt)*dt_sav(jat)

       dsi0d =   dpsi0(xd)*ddi_sav(iat)
       dsi1d =   dpsi1(xd)
       dsi2d =   dpsi2(xd)*dd_sav(iat)

       dsi0md = -dpsi0(mxd)*ddi_sav(iat)
       dsi1md =  dpsi1(mxd)
       dsi2md = -dpsi2(mxd)*dd_sav(iat)

       ! second derivatives of the weight functions
       ddsi0t =   ddpsi0(xt)*dt2i_sav(jat)
       ddsi1t =   ddpsi1(xt)*dti_sav(jat)
       ddsi2t =   ddpsi2(xt)

       ddsi0mt =  ddpsi0(mxt)*dt2i_sav(jat)
       ddsi1mt = -ddpsi1(mxt)*dti_sav(jat)
       ddsi2mt =  ddpsi2(mxt)

       !        ddsi0d =   ddpsi0(xd)*dd2i_sav(iat)
       !        ddsi1d =   ddpsi1(xd)*ddi_sav(iat)
       !        ddsi2d =   ddpsi2(xd)

       !        ddsi0md =  ddpsi0(mxd)*dd2i_sav(iat)
       !        ddsi1md = -ddpsi1(mxd)*ddi_sav(iat)
       !        ddsi2md =  ddpsi2(mxd)


       ! the free energy
       free  = h5(iat,jat, &
            si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
            si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

       ! derivative with respect to density
       df_d  = h5(iat,jat, &
            si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
            dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)


       ! derivative with respect to temperature
       df_t = h5(iat,jat, &
            dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
            si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

       ! derivative with respect to density**2
       !        df_dd = h5(iat,jat,
       !     1          si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
       !     2          ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md)

       ! derivative with respect to temperature**2
       df_tt = h5(iat,jat, &
            ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, &
            si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

       ! derivative with respect to temperature and density
       df_dt = h5(iat,jat, &
            dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
            dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)



       ! now get the pressure derivative with density, chemical potential, and
       ! electron positron number densities
       ! get the interpolation weight functions
       si0t   =  xpsi0(xt)
       si1t   =  xpsi1(xt)*dt_sav(jat)

       si0mt  =  xpsi0(mxt)
       si1mt  =  -xpsi1(mxt)*dt_sav(jat)

       si0d   =  xpsi0(xd)
       si1d   =  xpsi1(xd)*dd_sav(iat)

       si0md  =  xpsi0(mxd)
       si1md  =  -xpsi1(mxd)*dd_sav(iat)


       ! derivatives of weight functions
       dsi0t  = xdpsi0(xt)*dti_sav(jat)
       dsi1t  = xdpsi1(xt)

       dsi0mt = -xdpsi0(mxt)*dti_sav(jat)
       dsi1mt = xdpsi1(mxt)

       dsi0d  = xdpsi0(xd)*ddi_sav(iat)
       dsi1d  = xdpsi1(xd)

       dsi0md = -xdpsi0(mxd)*ddi_sav(iat)
       dsi1md = xdpsi1(mxd)


       ! look in the pressure derivative only once
       fi(1)  = dpdf(iat,jat)
       fi(2)  = dpdf(iat+1,jat)
       fi(3)  = dpdf(iat,jat+1)
       fi(4)  = dpdf(iat+1,jat+1)
       fi(5)  = dpdft(iat,jat)
       fi(6)  = dpdft(iat+1,jat)
       fi(7)  = dpdft(iat,jat+1)
       fi(8)  = dpdft(iat+1,jat+1)
       fi(9)  = dpdfd(iat,jat)
       fi(10) = dpdfd(iat+1,jat)
       fi(11) = dpdfd(iat,jat+1)
       fi(12) = dpdfd(iat+1,jat+1)
       fi(13) = dpdfdt(iat,jat)
       fi(14) = dpdfdt(iat+1,jat)
       fi(15) = dpdfdt(iat,jat+1)
       fi(16) = dpdfdt(iat+1,jat+1)

       ! pressure derivative with density
       dpepdd  = h3(iat,jat, &
            si0t,   si1t,   si0mt,   si1mt, &
            si0d,   si1d,   si0md,   si1md)
       dpepdd  = max(ye * dpepdd,1.0d-30)



       ! look in the electron chemical potential table only once
       fi(1)  = ef(iat,jat)
       fi(2)  = ef(iat+1,jat)
       fi(3)  = ef(iat,jat+1)
       fi(4)  = ef(iat+1,jat+1)
       fi(5)  = eft(iat,jat)
       fi(6)  = eft(iat+1,jat)
       fi(7)  = eft(iat,jat+1)
       fi(8)  = eft(iat+1,jat+1)
       fi(9)  = efd(iat,jat)
       fi(10) = efd(iat+1,jat)
       fi(11) = efd(iat,jat+1)
       fi(12) = efd(iat+1,jat+1)
       fi(13) = efdt(iat,jat)
       fi(14) = efdt(iat+1,jat)
       fi(15) = efdt(iat,jat+1)
       fi(16) = efdt(iat+1,jat+1)


       ! electron chemical potential etaele
       etaele  = h3(iat,jat, &
            si0t,   si1t,   si0mt,   si1mt, &
            si0d,   si1d,   si0md,   si1md)


       ! derivative with respect to density
       x       = h3(iat,jat, &
            si0t,   si1t,   si0mt,   si1mt, &
            dsi0d,  dsi1d,  dsi0md,  dsi1md)
       detadd  = ye * x

       ! derivative with respect to temperature
       detadt  = h3(iat,jat, &
            dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
            si0d,   si1d,   si0md,   si1md)

       ! derivative with respect to abar and zbar
       detada = -x * din * ytot1
       detadz =  x * den * ytot1



       ! look in the number density table only once
       fi(1)  = xf(iat,jat)
       fi(2)  = xf(iat+1,jat)
       fi(3)  = xf(iat,jat+1)
       fi(4)  = xf(iat+1,jat+1)
       fi(5)  = xft(iat,jat)
       fi(6)  = xft(iat+1,jat)
       fi(7)  = xft(iat,jat+1)
       fi(8)  = xft(iat+1,jat+1)
       fi(9)  = xfd(iat,jat)
       fi(10) = xfd(iat+1,jat)
       fi(11) = xfd(iat,jat+1)
       fi(12) = xfd(iat+1,jat+1)
       fi(13) = xfdt(iat,jat)
       fi(14) = xfdt(iat+1,jat)
       fi(15) = xfdt(iat,jat+1)
       fi(16) = xfdt(iat+1,jat+1)

       ! electron + positron number densities
       xnefer   = h3(iat,jat, &
            si0t,   si1t,   si0mt,   si1mt, &
            si0d,   si1d,   si0md,   si1md)

       ! derivative with respect to density
       x        = h3(iat,jat, &
            si0t,   si1t,   si0mt,   si1mt, &
            dsi0d,  dsi1d,  dsi0md,  dsi1md)
       x = max(x,1.0d-30)
       dxnedd   = ye * x

       ! derivative with respect to temperature
       dxnedt   = h3(iat,jat, &
            dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
            si0d,   si1d,   si0md,   si1md)

       ! derivative with respect to abar and zbar
       dxneda = -x * din * ytot1
       dxnedz =  x  * den * ytot1


       ! the desired electron-positron thermodynamic quantities

       ! dpepdd at high temperatures and low densities is below the
       ! floating point limit of the subtraction of two large terms.
       ! since dpresdd doesn't enter the maxwell relations at all, use the
       ! bicubic interpolation done above instead of the formally correct expression
       x       = din * din
       pele    = x * df_d
       dpepdt  = x * df_dt
       !        dpepdd  = ye * (x * df_dd + 2.0d0 * din * df_d)
       s       = dpepdd/ye - 2.0d0 * din * df_d
       dpepda  = -ytot1 * (2.0d0 * pele + s * din)
       dpepdz  = den*ytot1*(2.0d0 * din * df_d  +  s)


       x       = ye * ye
       sele    = -df_t * ye
       dsepdt  = -df_tt * ye
       dsepdd  = -df_dt * x
       dsepda  = ytot1 * (ye * df_dt * din - sele)
       dsepdz  = -ytot1 * (ye * df_dt * den  + df_t)


       eele    = ye*free + temp * sele
       deepdt  = temp * dsepdt
       deepdd  = x * df_d + temp * dsepdd
       deepda  = -ye * ytot1 * (free +  df_d * din) + temp * dsepda
       deepdz  = ytot1* (free + ye * df_d * den) + temp * dsepdz




       ! coulomb section:

       ! uniform background corrections only
       ! from yakovlev & shalybkov 1989
       ! lami is the average ion seperation
       ! plasg is the plasma coupling parameter

       z        = forth * pi
       s        = z * xni
       dsdd     = z * dxnidd
       dsda     = z * dxnida

       lami     = 1.0d0/s**third
       inv_lami = 1.0d0/lami
       z        = -third * lami
       lamidd   = z * dsdd/s
       lamida   = z * dsda/s

       plasg    = zbar*zbar*esqu*ktinv*inv_lami
       z        = -plasg * inv_lami
       plasgdd  = z * lamidd
       plasgda  = z * lamida
       plasgdt  = -plasg*ktinv * kerg
       plasgdz  = 2.0d0 * plasg/zbar


       ! yakovlev & shalybkov 1989 equations 82, 85, 86, 87
       if (plasg .ge. 1.0) then
          x        = plasg**(0.25d0)
          y        = avo * ytot1 * kerg
          ecoul    = y * temp * (a1*plasg + b1*x + c1/x + d1)
          pcoul    = third * den * ecoul
          scoul    = -y * (3.0d0*b1*x - 5.0d0*c1/x &
               + d1 * (log(plasg) - 1.0d0) - e1)

          y        = avo*ytot1*kt*(a1 + 0.25d0/plasg*(b1*x - c1/x))
          decouldd = y * plasgdd
          decouldt = y * plasgdt + ecoul/temp
          decoulda = y * plasgda - ecoul/abar
          decouldz = y * plasgdz

          y        = third * den
          dpcouldd = third * ecoul + y*decouldd
          dpcouldt = y * decouldt
          dpcoulda = y * decoulda
          dpcouldz = y * decouldz


          y        = -avo*kerg/(abar*plasg)*(0.75d0*b1*x+1.25d0*c1/x+d1)
          dscouldd = y * plasgdd
          dscouldt = y * plasgdt
          dscoulda = y * plasgda - scoul/abar
          dscouldz = y * plasgdz


          ! yakovlev & shalybkov 1989 equations 102, 103, 104
       else if (plasg .lt. 1.0) then
          x        = plasg*sqrt(plasg)
          y        = plasg**b2
          z        = c2 * x - third * a2 * y
          pcoul    = -pion * z
          ecoul    = 3.0d0 * pcoul/den
          scoul    = -avo/abar*kerg*(c2*x -a2*(b2-1.0d0)/b2*y)

          s        = 1.5d0*c2*x/plasg - third*a2*b2*y/plasg
          dpcouldd = -dpiondd*z - pion*s*plasgdd
          dpcouldt = -dpiondt*z - pion*s*plasgdt
          dpcoulda = -dpionda*z - pion*s*plasgda
          dpcouldz = -dpiondz*z - pion*s*plasgdz

          s        = 3.0d0/den
          decouldd = s * dpcouldd - ecoul/den
          decouldt = s * dpcouldt
          decoulda = s * dpcoulda
          decouldz = s * dpcouldz

          s        = -avo*kerg/(abar*plasg)*(1.5d0*c2*x-a2*(b2-1.0d0)*y)
          dscouldd = s * plasgdd
          dscouldt = s * plasgdt
          dscoulda = s * plasgda - scoul/abar
          dscouldz = s * plasgdz
       end if


       ! bomb proof
       x   = prad + pion + pele + pcoul
       y   = erad + eion + eele + ecoul
       z   = srad + sion + sele + scoul

       !        write(6,*) x,y,z
       !        if (x .le. 0.0 .or. y .le. 0.0 .or. z .le. 0.0) then
       if (x .le. 0.0 .or. y .le. 0.0) then
          !        if (x .le. 0.0) then

          !         write(6,*)
          !         write(6,*) 'coulomb corrections are causing a negative pressure'
          !         write(6,*) 'setting all coulomb corrections to zero'
          !         write(6,*)

          pcoul    = 0.0d0
          dpcouldd = 0.0d0
          dpcouldt = 0.0d0
          dpcoulda = 0.0d0
          dpcouldz = 0.0d0
          ecoul    = 0.0d0
          decouldd = 0.0d0
          decouldt = 0.0d0
          decoulda = 0.0d0
          decouldz = 0.0d0
          scoul    = 0.0d0
          dscouldd = 0.0d0
          dscouldt = 0.0d0
          dscoulda = 0.0d0
          dscouldz = 0.0d0
       end if


       ! sum all the gas components
       pgas    = pion + pele + pcoul
       egas    = eion + eele + ecoul
       sgas    = sion + sele + scoul

       dpgasdd = dpiondd + dpepdd + dpcouldd
       dpgasdt = dpiondt + dpepdt + dpcouldt
       dpgasda = dpionda + dpepda + dpcoulda
       dpgasdz = dpiondz + dpepdz + dpcouldz

       degasdd = deiondd + deepdd + decouldd
       degasdt = deiondt + deepdt + decouldt
       degasda = deionda + deepda + decoulda
       degasdz = deiondz + deepdz + decouldz

       dsgasdd = dsiondd + dsepdd + dscouldd
       dsgasdt = dsiondt + dsepdt + dscouldt
       dsgasda = dsionda + dsepda + dscoulda
       dsgasdz = dsiondz + dsepdz + dscouldz




       ! add in radiation to get the total
       pres    = prad + pgas
       ener    = erad + egas
       entr    = srad + sgas

       dpresdd = dpraddd + dpgasdd
       dpresdt = dpraddt + dpgasdt
       dpresda = dpradda + dpgasda
       dpresdz = dpraddz + dpgasdz

       denerdd = deraddd + degasdd
       denerdt = deraddt + degasdt
       denerda = deradda + degasda
       denerdz = deraddz + degasdz

       dentrdd = dsraddd + dsgasdd
       dentrdt = dsraddt + dsgasdt
       dentrda = dsradda + dsgasda
       dentrdz = dsraddz + dsgasdz


       ! for the gas
       ! the temperature and density exponents (c&g 9.81 9.82)
       ! the specific heat at constant volume (c&g 9.92)
       ! the third adiabatic exponent (c&g 9.93)
       ! the first adiabatic exponent (c&g 9.97)
       ! the second adiabatic exponent (c&g 9.105)
       ! the specific heat at constant pressure (c&g 9.98)
       ! and relativistic formula for the sound speed (c&g 14.29)

       zz        = pgas*deni
       zzi       = den/pgas
       chit_gas  = temp/pgas * dpgasdt
       chid_gas  = dpgasdd*zzi
       cv_gas    = degasdt
       x         = zz * chit_gas/(temp * cv_gas)
       gam3_gas  = x + 1.0d0
       gam1_gas  = chit_gas*x + chid_gas
       nabad_gas = x/gam1_gas
       gam2_gas  = 1.0d0/(1.0d0 - nabad_gas)
       cp_gas    = cv_gas * gam1_gas/chid_gas
       z         = 1.0d0 + (egas + light2)*zzi
       sound_gas = clight * sqrt(gam1_gas/z)



       ! for the totals
       zz    = pres*deni
       zzi   = den/pres
       chit  = temp/pres * dpresdt
       chid  = dpresdd*zzi
       cv    = denerdt
       x     = zz * chit/(temp * cv)
       gam3  = x + 1.0d0
       gam1  = chit*x + chid
       nabad = x/gam1
       gam2  = 1.0d0/(1.0d0 - nabad)
       cp    = cv * gam1/chid
       z     = 1.0d0 + (ener + light2)*zzi
       sound = clight * sqrt(gam1/z)



       ! maxwell relations; each is zero if the consistency is perfect
       x   = den * den

       dse = temp*dentrdt/denerdt - 1.0d0

       dpe = (denerdd*x + temp*dpresdt)/pres - 1.0d0

       dsp = -dentrdd*x/dpresdt - 1.0d0


       ! store this row
       ptot_row(j)   = pres
       dpt_row(j)    = dpresdt
       dpd_row(j)    = dpresdd
       dpa_row(j)    = dpresda
       dpz_row(j)    = dpresdz

       etot_row(j)   = ener
       det_row(j)    = denerdt
       ded_row(j)    = denerdd
       dea_row(j)    = denerda
       dez_row(j)    = denerdz

       stot_row(j)   = entr
       dst_row(j)    = dentrdt
       dsd_row(j)    = dentrdd
       dsa_row(j)    = dentrda
       dsz_row(j)    = dentrdz


       pgas_row(j)   = pgas
       dpgast_row(j) = dpgasdt
       dpgasd_row(j) = dpgasdd
       dpgasa_row(j) = dpgasda
       dpgasz_row(j) = dpgasdz

       egas_row(j)   = egas
       degast_row(j) = degasdt
       degasd_row(j) = degasdd
       degasa_row(j) = degasda
       degasz_row(j) = degasdz

       sgas_row(j)   = sgas
       dsgast_row(j) = dsgasdt
       dsgasd_row(j) = dsgasdd
       dsgasa_row(j) = dsgasda
       dsgasz_row(j) = dsgasdz


       prad_row(j)   = prad
       dpradt_row(j) = dpraddt
       dpradd_row(j) = dpraddd
       dprada_row(j) = dpradda
       dpradz_row(j) = dpraddz

       erad_row(j)   = erad
       deradt_row(j) = deraddt
       deradd_row(j) = deraddd
       derada_row(j) = deradda
       deradz_row(j) = deraddz

       srad_row(j)   = srad
       dsradt_row(j) = dsraddt
       dsradd_row(j) = dsraddd
       dsrada_row(j) = dsradda
       dsradz_row(j) = dsraddz


       pion_row(j)   = pion
       dpiont_row(j) = dpiondt
       dpiond_row(j) = dpiondd
       dpiona_row(j) = dpionda
       dpionz_row(j) = dpiondz

       eion_row(j)   = eion
       deiont_row(j) = deiondt
       deiond_row(j) = deiondd
       deiona_row(j) = deionda
       deionz_row(j) = deiondz

       sion_row(j)   = sion
       dsiont_row(j) = dsiondt
       dsiond_row(j) = dsiondd
       dsiona_row(j) = dsionda
       dsionz_row(j) = dsiondz

       xni_row(j)    = xni

       pele_row(j)   = pele
       ppos_row(j)   = 0.0d0
       dpept_row(j)  = dpepdt
       dpepd_row(j)  = dpepdd
       dpepa_row(j)  = dpepda
       dpepz_row(j)  = dpepdz

       eele_row(j)   = eele
       epos_row(j)   = 0.0d0
       deept_row(j)  = deepdt
       deepd_row(j)  = deepdd
       deepa_row(j)  = deepda
       deepz_row(j)  = deepdz

       sele_row(j)   = sele
       spos_row(j)   = 0.0d0
       dsept_row(j)  = dsepdt
       dsepd_row(j)  = dsepdd
       dsepa_row(j)  = dsepda
       dsepz_row(j)  = dsepdz

       xnem_row(j)   = xnem
       xne_row(j)    = xnefer
       dxnet_row(j)  = dxnedt
       dxned_row(j)  = dxnedd
       dxnea_row(j)  = dxneda
       dxnez_row(j)  = dxnedz
       xnp_row(j)    = 0.0d0
       zeff_row(j)   = zbar

       etaele_row(j) = etaele
       detat_row(j)  = detadt
       detad_row(j)  = detadd
       detaa_row(j)  = detada
       detaz_row(j)  = detadz
       etapos_row(j) = 0.0d0

       pcou_row(j)   = pcoul
       dpcout_row(j) = dpcouldt
       dpcoud_row(j) = dpcouldd
       dpcoua_row(j) = dpcoulda
       dpcouz_row(j) = dpcouldz

       ecou_row(j)   = ecoul
       decout_row(j) = decouldt
       decoud_row(j) = decouldd
       decoua_row(j) = decoulda
       decouz_row(j) = decouldz

       scou_row(j)   = scoul
       dscout_row(j) = dscouldt
       dscoud_row(j) = dscouldd
       dscoua_row(j) = dscoulda
       dscouz_row(j) = dscouldz

       plasg_row(j)  = plasg

       dse_row(j)    = dse
       dpe_row(j)    = dpe
       dsp_row(j)    = dsp

       cv_gas_row(j)    = cv_gas
       cp_gas_row(j)    = cp_gas
       gam1_gas_row(j)  = gam1_gas
       gam2_gas_row(j)  = gam2_gas
       gam3_gas_row(j)  = gam3_gas
       nabad_gas_row(j) = nabad_gas
       cs_gas_row(j)    = sound_gas

       cv_row(j)     = cv
       cp_row(j)     = cp
       gam1_row(j)   = gam1
       gam2_row(j)   = gam2
       gam3_row(j)   = gam3
       nabad_row(j)  = nabad
       cs_row(j)     = sound

       ! end of pipeline loop
    enddo
    return
  end subroutine helmeos

  subroutine eos_get_misc(rho_in, temp_in, ye_in, ytot_in, mexc_in, &
       eta_out)
    real(8),intent(in) :: rho_in, temp_in, ye_in, ytot_in, mexc_in
    real(8),intent(out) :: eta_out

    

    ! declare
    integer          i,j
    double precision temp,den,abar,zbar,ytot1,ye, &
         x,y,zz,zzi,deni,tempi,xni,dxnidd,dxnida, &
         dpepdt,dpepdd,deepdt,deepdd,dsepdd,dsepdt, &
         dpraddd,dpraddt,deraddd,deraddt,dpiondd,dpiondt, &
         deiondd,deiondt,dsraddd,dsraddt,dsiondd,dsiondt, &
         dse,dpe,dsp,kt,ktinv,prad,erad,srad,pion,eion, &
         sion,xnem,pele,eele,sele,pres,ener,entr,dpresdd, &
         dpresdt,denerdd,denerdt,dentrdd,dentrdt,cv,cp, &
         gam1,gam2,gam3,chit,chid,nabad,sound,etaele, &
         detadt,detadd,xnefer,dxnedt,dxnedd,s

    double precision pgas,dpgasdd,dpgasdt,dpgasda,dpgasdz, &
         egas,degasdd,degasdt,degasda,degasdz, &
         sgas,dsgasdd,dsgasdt,dsgasda,dsgasdz, &
         cv_gas,cp_gas,gam1_gas,gam2_gas,gam3_gas, &
         chit_gas,chid_gas,nabad_gas,sound_gas


    double precision sioncon,forth,forpi,kergavo,ikavo,asoli3,light2
    parameter        (sioncon = (2.0d0 * pi * amu * kerg)/(h*h), &
         forth   = 4.0d0/3.0d0, &
         forpi   = 4.0d0 * pi, &
         kergavo = kerg * avo, &
         ikavo   = 1.0d0/kergavo, &
         asoli3  = asol/3.0d0, &
         light2  = clight * clight)

    ! for the abar derivatives
    double precision dpradda,deradda,dsradda, &
         dpionda,deionda,dsionda, &
         dpepda,deepda,dsepda, &
         dpresda,denerda,dentrda, &
         detada,dxneda

    ! for the zbar derivatives
    double precision dpraddz,deraddz,dsraddz, &
         dpiondz,deiondz,dsiondz, &
         dpepdz,deepdz,dsepdz, &
         dpresdz,denerdz,dentrdz, &
         detadz,dxnedz

    ! for the interpolations
    integer          iat,jat
    double precision free,df_d,df_t,df_dd,df_tt,df_dt
    double precision xt,xd,mxt,mxd, &
         si0t,si1t,si2t,si0mt,si1mt,si2mt, &
         si0d,si1d,si2d,si0md,si1md,si2md, &
         dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt, &
         dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md, &
         ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt, &
         ddsi0d,ddsi1d,ddsi2d,ddsi0md,ddsi1md,ddsi2md, &
         z,psi0,dpsi0,ddpsi0,psi1,dpsi1,ddpsi1,psi2, &
         dpsi2,ddpsi2,din,h5,fi(36), &
         xpsi0,xdpsi0,xpsi1,xdpsi1,h3, &
         w0t,w1t,w2t,w0mt,w1mt,w2mt, &
         w0d,w1d,w2d,w0md,w1md,w2md


    ! for the uniform background coulomb correction
    double precision dsdd,dsda,lami,inv_lami,lamida,lamidd, &
         plasg,plasgdd,plasgdt,plasgda,plasgdz, &
         ecoul,decouldd,decouldt,decoulda,decouldz, &
         pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz, &
         scoul,dscouldd,dscouldt,dscoulda,dscouldz, &
         a1,b1,c1,d1,e1,a2,b2,c2,third,esqu
    parameter        (a1    = -0.898004d0, &
         b1    =  0.96786d0, &
         c1    =  0.220703d0, &
         d1    = -0.86097d0, &
         e1    =  2.5269d0, &
         a2    =  0.29561d0, &
         b2    =  1.9885d0, &
         c2    =  0.288675d0, &
         third =  1.0d0/3.0d0, &
         esqu  =  qe * qe)


    ! quintic hermite polynomial statement functions
    ! psi0 and its derivatives
    psi0(z)   = z**3 * ( z * (-6.0d0*z + 15.0d0) -10.0d0) + 1.0d0
    dpsi0(z)  = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
    ddpsi0(z) = z* ( z*( -120.0d0*z + 180.0d0) -60.0d0)


    ! psi1 and its derivatives
    psi1(z)   = z* ( z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
    dpsi1(z)  = z*z * ( z * (-15.0d0*z + 32.0d0) - 18.0d0) +1.0d0
    ddpsi1(z) = z * (z * (-60.0d0*z + 96.0d0) -36.0d0)


    ! psi2  and its derivatives
    psi2(z)   = 0.5d0*z*z*( z* ( z * (-z + 3.0d0) - 3.0d0) + 1.0d0)
    dpsi2(z)  = 0.5d0*z*( z*(z*(-5.0d0*z + 12.0d0) - 9.0d0) + 2.0d0)
    ddpsi2(z) = 0.5d0*(z*( z * (-20.0d0*z + 36.0d0) - 18.0d0) + 2.0d0)


    ! biquintic hermite polynomial statement function
    h5(i,j,w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md)= &
         fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t &
         + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt &
         + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t &
         + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt &
         + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t &
         + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt &
         + fi(13) *w1d*w0t   + fi(14) *w1md*w0t &
         + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt &
         + fi(17) *w2d*w0t   + fi(18) *w2md*w0t &
         + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt &
         + fi(21) *w1d*w1t   + fi(22) *w1md*w1t &
         + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt &
         + fi(25) *w2d*w1t   + fi(26) *w2md*w1t &
         + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt &
         + fi(29) *w1d*w2t   + fi(30) *w1md*w2t &
         + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt &
         + fi(33) *w2d*w2t   + fi(34) *w2md*w2t &
         + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt



    ! cubic hermite polynomial statement functions
    ! psi0 & derivatives
    xpsi0(z)  = z * z * (2.0d0*z - 3.0d0) + 1.0
    xdpsi0(z) = z * (6.0d0*z - 6.0d0)


    ! psi1 & derivatives
    xpsi1(z)  = z * ( z * (z - 2.0d0) + 1.0d0)
    xdpsi1(z) = z * (3.0d0*z - 4.0d0) + 1.0d0


    ! bicubic hermite polynomial statement function
    h3(i,j,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) = &
         fi(1)  *w0d*w0t   +  fi(2)  *w0md*w0t &
         + fi(3)  *w0d*w0mt  +  fi(4)  *w0md*w0mt &
         + fi(5)  *w0d*w1t   +  fi(6)  *w0md*w1t &
         + fi(7)  *w0d*w1mt  +  fi(8)  *w0md*w1mt &
         + fi(9)  *w1d*w0t   +  fi(10) *w1md*w0t &
         + fi(11) *w1d*w0mt  +  fi(12) *w1md*w0mt &
         + fi(13) *w1d*w1t   +  fi(14) *w1md*w1t &
         + fi(15) *w1d*w1mt  +  fi(16) *w1md*w1mt



    ! popular format statements
01  format(1x,5(a,1pe11.3))
02  format(1x,a,1p4e16.8)
03  format(1x,4(a,1pe11.3))
04  format(1x,4(a,i4))

    ! start of pipeline loop, normal execution starts here
    ! eosfail = .false.

    !       if (temp_row(j) .le. 0.0) stop 'temp less than 0 in helmeos'
    !       if (den_row(j)  .le. 0.0) stop 'den less than 0 in helmeos'

    temp  = temp_in
    den   = rho_in
    abar  = 1d0/ytot_in
    zbar  = abar*ye_in
    ytot1 = 1.0d0/abar
    ye    = ye_in

    ! initialize
    deni    = 1.0d0/den
    tempi   = 1.0d0/temp
    kt      = kerg * temp
    ktinv   = 1.0d0/kt


    ! radiation section:
    prad    = asoli3 * temp * temp * temp * temp
    dpraddd = 0.0d0
    dpraddt = 4.0d0 * prad*tempi
    dpradda = 0.0d0
    dpraddz = 0.0d0

    erad    = 3.0d0 * prad*deni
    deraddd = -erad*deni
    deraddt = 3.0d0 * dpraddt*deni
    deradda = 0.0d0
    deraddz = 0.0d0

    srad    = (prad*deni + erad)*tempi
    dsraddd = (dpraddd*deni - prad*deni*deni + deraddd)*tempi
    dsraddt = (dpraddt*deni + deraddt - srad)*tempi
    dsradda = 0.0d0
    dsraddz = 0.0d0


    ! ion section:
    xni     = avo * ytot1 * den
    dxnidd  = avo * ytot1
    dxnida  = -xni * ytot1

    pion    = xni * kt
    dpiondd = dxnidd * kt
    dpiondt = xni * kerg
    dpionda = dxnida * kt
    dpiondz = 0.0d0

    eion    = 1.5d0 * pion*deni
    deiondd = (1.5d0 * dpiondd - eion)*deni
    deiondt = 1.5d0 * dpiondt*deni
    deionda = 1.5d0 * dpionda*deni
    deiondz = 0.0d0


    ! sackur-tetrode equation for the ion entropy of
    ! a single ideal gas characterized by abar
    x       = abar*abar*sqrt(abar) * deni/avo
    s       = sioncon * temp
    z       = x * s * sqrt(s)
    y       = log(z)

    !        y       = 1.0d0/(abar*kt)
    !        yy      = y * sqrt(y)
    !        z       = xni * sifac * yy
    !        etaion  = log(z)


    sion    = (pion*deni + eion)*tempi + kergavo * ytot1 * y
    dsiondd = (dpiondd*deni - pion*deni*deni + deiondd)*tempi &
         - kergavo * deni * ytot1
    dsiondt = (dpiondt*deni + deiondt)*tempi - &
         (pion*deni + eion) * tempi*tempi &
         + 1.5d0 * kergavo * tempi*ytot1
    x       = avo*kerg/abar
    dsionda = (dpionda*deni + deionda)*tempi &
         + kergavo*ytot1*ytot1* (2.5d0 - y)
    dsiondz = 0.0d0



    ! electron-positron section:


    ! assume complete ionization
    xnem    = xni * zbar


    ! enter the table with ye*den
    din = ye*den


    ! bomb proof the input
    temp = min( temp, t(jmax))
    temp = max( temp, t(1))
    din  = min( din , d(imax))
    din  = max( din , d(1))

    ! if (temp .gt. t(jmax)) then
    !    write(6,01) 'temp=',temp,' t(jmax)=',t(jmax)
    !    write(6,*) 'temp too hot, off grid'
    !    write(6,*) 'setting eosfail to true and returning'
    !    eosfail = .true.
    !    return
    ! end if
    ! if (temp .lt. t(1)) then
    !    write(6,01) 'temp=',temp,' t(1)=',t(1)
    !    write(6,*) 'temp too cold, off grid'
    !    write(6,*) 'setting eosfail to true and returning'
    !    eosfail = .true.
    !    return
    ! end if
    ! if (din  .gt. d(imax)) then
    !    write(6,01) 'den*ye=',din,' d(imax)=',d(imax)
    !    write(6,*) 'ye*den too big, off grid'
    !    write(6,*) 'setting eosfail to true and returning'
    !    eosfail = .true.
    !    return
    ! end if
    ! if (din  .lt. d(1)) then
    !    write(6,01) 'ye*den=',din,' d(1)=',d(1)
    !    write(6,*) 'ye*den too small, off grid'
    !    write(6,*) 'setting eosfail to true and returning'
    !    eosfail = .true.
    !    return
    ! end if

    ! hash locate this temperature and density
    jat = int((log10(temp) - tlo)*tstpi) + 1
    jat = max(1,min(jat,jmax-1))
    iat = int((log10(din) - dlo)*dstpi) + 1
    iat = max(1,min(iat,imax-1))


    ! access the table locations only once
    fi(1)  = f(iat,jat)
    fi(2)  = f(iat+1,jat)
    fi(3)  = f(iat,jat+1)
    fi(4)  = f(iat+1,jat+1)
    fi(5)  = ft(iat,jat)
    fi(6)  = ft(iat+1,jat)
    fi(7)  = ft(iat,jat+1)
    fi(8)  = ft(iat+1,jat+1)
    fi(9)  = ftt(iat,jat)
    fi(10) = ftt(iat+1,jat)
    fi(11) = ftt(iat,jat+1)
    fi(12) = ftt(iat+1,jat+1)
    fi(13) = fd(iat,jat)
    fi(14) = fd(iat+1,jat)
    fi(15) = fd(iat,jat+1)
    fi(16) = fd(iat+1,jat+1)
    fi(17) = fdd(iat,jat)
    fi(18) = fdd(iat+1,jat)
    fi(19) = fdd(iat,jat+1)
    fi(20) = fdd(iat+1,jat+1)
    fi(21) = fdt(iat,jat)
    fi(22) = fdt(iat+1,jat)
    fi(23) = fdt(iat,jat+1)
    fi(24) = fdt(iat+1,jat+1)
    fi(25) = fddt(iat,jat)
    fi(26) = fddt(iat+1,jat)
    fi(27) = fddt(iat,jat+1)
    fi(28) = fddt(iat+1,jat+1)
    fi(29) = fdtt(iat,jat)
    fi(30) = fdtt(iat+1,jat)
    fi(31) = fdtt(iat,jat+1)
    fi(32) = fdtt(iat+1,jat+1)
    fi(33) = fddtt(iat,jat)
    fi(34) = fddtt(iat+1,jat)
    fi(35) = fddtt(iat,jat+1)
    fi(36) = fddtt(iat+1,jat+1)


    ! various differences
    xt  = max( (temp - t(jat))*dti_sav(jat), 0.0d0)
    xd  = max( (din - d(iat))*ddi_sav(iat), 0.0d0)
    mxt = 1.0d0 - xt
    mxd = 1.0d0 - xd

    ! the six density and six temperature basis functions
    si0t =   psi0(xt)
    si1t =   psi1(xt)*dt_sav(jat)
    si2t =   psi2(xt)*dt2_sav(jat)

    si0mt =  psi0(mxt)
    si1mt = -psi1(mxt)*dt_sav(jat)
    si2mt =  psi2(mxt)*dt2_sav(jat)

    si0d =   psi0(xd)
    si1d =   psi1(xd)*dd_sav(iat)
    si2d =   psi2(xd)*dd2_sav(iat)

    si0md =  psi0(mxd)
    si1md = -psi1(mxd)*dd_sav(iat)
    si2md =  psi2(mxd)*dd2_sav(iat)

    ! derivatives of the weight functions
    dsi0t =   dpsi0(xt)*dti_sav(jat)
    dsi1t =   dpsi1(xt)
    dsi2t =   dpsi2(xt)*dt_sav(jat)

    dsi0mt = -dpsi0(mxt)*dti_sav(jat)
    dsi1mt =  dpsi1(mxt)
    dsi2mt = -dpsi2(mxt)*dt_sav(jat)

    dsi0d =   dpsi0(xd)*ddi_sav(iat)
    dsi1d =   dpsi1(xd)
    dsi2d =   dpsi2(xd)*dd_sav(iat)

    dsi0md = -dpsi0(mxd)*ddi_sav(iat)
    dsi1md =  dpsi1(mxd)
    dsi2md = -dpsi2(mxd)*dd_sav(iat)

    ! second derivatives of the weight functions
    ddsi0t =   ddpsi0(xt)*dt2i_sav(jat)
    ddsi1t =   ddpsi1(xt)*dti_sav(jat)
    ddsi2t =   ddpsi2(xt)

    ddsi0mt =  ddpsi0(mxt)*dt2i_sav(jat)
    ddsi1mt = -ddpsi1(mxt)*dti_sav(jat)
    ddsi2mt =  ddpsi2(mxt)

    !        ddsi0d =   ddpsi0(xd)*dd2i_sav(iat)
    !        ddsi1d =   ddpsi1(xd)*ddi_sav(iat)
    !        ddsi2d =   ddpsi2(xd)

    !        ddsi0md =  ddpsi0(mxd)*dd2i_sav(iat)
    !        ddsi1md = -ddpsi1(mxd)*ddi_sav(iat)
    !        ddsi2md =  ddpsi2(mxd)


    ! the free energy
    free  = h5(iat,jat, &
         si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
         si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

    ! derivative with respect to density
    df_d  = h5(iat,jat, &
         si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
         dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)


    ! derivative with respect to temperature
    df_t = h5(iat,jat, &
         dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
         si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

    ! derivative with respect to density**2
    !        df_dd = h5(iat,jat,
    !     1          si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
    !     2          ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md)

    ! derivative with respect to temperature**2
    df_tt = h5(iat,jat, &
         ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, &
         si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

    ! derivative with respect to temperature and density
    df_dt = h5(iat,jat, &
         dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
         dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)



    ! now get the pressure derivative with density, chemical potential, and
    ! electron positron number densities
    ! get the interpolation weight functions
    si0t   =  xpsi0(xt)
    si1t   =  xpsi1(xt)*dt_sav(jat)

    si0mt  =  xpsi0(mxt)
    si1mt  =  -xpsi1(mxt)*dt_sav(jat)

    si0d   =  xpsi0(xd)
    si1d   =  xpsi1(xd)*dd_sav(iat)

    si0md  =  xpsi0(mxd)
    si1md  =  -xpsi1(mxd)*dd_sav(iat)


    ! derivatives of weight functions
    dsi0t  = xdpsi0(xt)*dti_sav(jat)
    dsi1t  = xdpsi1(xt)

    dsi0mt = -xdpsi0(mxt)*dti_sav(jat)
    dsi1mt = xdpsi1(mxt)

    dsi0d  = xdpsi0(xd)*ddi_sav(iat)
    dsi1d  = xdpsi1(xd)

    dsi0md = -xdpsi0(mxd)*ddi_sav(iat)
    dsi1md = xdpsi1(mxd)


    ! look in the pressure derivative only once
    fi(1)  = dpdf(iat,jat)
    fi(2)  = dpdf(iat+1,jat)
    fi(3)  = dpdf(iat,jat+1)
    fi(4)  = dpdf(iat+1,jat+1)
    fi(5)  = dpdft(iat,jat)
    fi(6)  = dpdft(iat+1,jat)
    fi(7)  = dpdft(iat,jat+1)
    fi(8)  = dpdft(iat+1,jat+1)
    fi(9)  = dpdfd(iat,jat)
    fi(10) = dpdfd(iat+1,jat)
    fi(11) = dpdfd(iat,jat+1)
    fi(12) = dpdfd(iat+1,jat+1)
    fi(13) = dpdfdt(iat,jat)
    fi(14) = dpdfdt(iat+1,jat)
    fi(15) = dpdfdt(iat,jat+1)
    fi(16) = dpdfdt(iat+1,jat+1)

    ! pressure derivative with density
    dpepdd  = h3(iat,jat, &
         si0t,   si1t,   si0mt,   si1mt, &
         si0d,   si1d,   si0md,   si1md)
    dpepdd  = max(ye * dpepdd,1.0d-30)



    ! look in the electron chemical potential table only once
    fi(1)  = ef(iat,jat)
    fi(2)  = ef(iat+1,jat)
    fi(3)  = ef(iat,jat+1)
    fi(4)  = ef(iat+1,jat+1)
    fi(5)  = eft(iat,jat)
    fi(6)  = eft(iat+1,jat)
    fi(7)  = eft(iat,jat+1)
    fi(8)  = eft(iat+1,jat+1)
    fi(9)  = efd(iat,jat)
    fi(10) = efd(iat+1,jat)
    fi(11) = efd(iat,jat+1)
    fi(12) = efd(iat+1,jat+1)
    fi(13) = efdt(iat,jat)
    fi(14) = efdt(iat+1,jat)
    fi(15) = efdt(iat,jat+1)
    fi(16) = efdt(iat+1,jat+1)


    ! electron chemical potential etaele
    etaele  = h3(iat,jat, &
         si0t,   si1t,   si0mt,   si1mt, &
         si0d,   si1d,   si0md,   si1md)


    ! derivative with respect to density
    x       = h3(iat,jat, &
         si0t,   si1t,   si0mt,   si1mt, &
         dsi0d,  dsi1d,  dsi0md,  dsi1md)
    detadd  = ye * x

    ! derivative with respect to temperature
    detadt  = h3(iat,jat, &
         dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
         si0d,   si1d,   si0md,   si1md)

    ! derivative with respect to abar and zbar
    detada = -x * din * ytot1
    detadz =  x * den * ytot1



    ! look in the number density table only once
    fi(1)  = xf(iat,jat)
    fi(2)  = xf(iat+1,jat)
    fi(3)  = xf(iat,jat+1)
    fi(4)  = xf(iat+1,jat+1)
    fi(5)  = xft(iat,jat)
    fi(6)  = xft(iat+1,jat)
    fi(7)  = xft(iat,jat+1)
    fi(8)  = xft(iat+1,jat+1)
    fi(9)  = xfd(iat,jat)
    fi(10) = xfd(iat+1,jat)
    fi(11) = xfd(iat,jat+1)
    fi(12) = xfd(iat+1,jat+1)
    fi(13) = xfdt(iat,jat)
    fi(14) = xfdt(iat+1,jat)
    fi(15) = xfdt(iat,jat+1)
    fi(16) = xfdt(iat+1,jat+1)

    ! electron + positron number densities
    xnefer   = h3(iat,jat, &
         si0t,   si1t,   si0mt,   si1mt, &
         si0d,   si1d,   si0md,   si1md)

    ! derivative with respect to density
    x        = h3(iat,jat, &
         si0t,   si1t,   si0mt,   si1mt, &
         dsi0d,  dsi1d,  dsi0md,  dsi1md)
    x = max(x,1.0d-30)
    dxnedd   = ye * x

    ! derivative with respect to temperature
    dxnedt   = h3(iat,jat, &
         dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
         si0d,   si1d,   si0md,   si1md)

    ! derivative with respect to abar and zbar
    dxneda = -x * din * ytot1
    dxnedz =  x  * den * ytot1


    ! the desired electron-positron thermodynamic quantities

    ! dpepdd at high temperatures and low densities is below the
    ! floating point limit of the subtraction of two large terms.
    ! since dpresdd doesn't enter the maxwell relations at all, use the
    ! bicubic interpolation done above instead of the formally correct expression
    x       = din * din
    pele    = x * df_d
    dpepdt  = x * df_dt
    !        dpepdd  = ye * (x * df_dd + 2.0d0 * din * df_d)
    s       = dpepdd/ye - 2.0d0 * din * df_d
    dpepda  = -ytot1 * (2.0d0 * pele + s * din)
    dpepdz  = den*ytot1*(2.0d0 * din * df_d  +  s)


    x       = ye * ye
    sele    = -df_t * ye
    dsepdt  = -df_tt * ye
    dsepdd  = -df_dt * x
    dsepda  = ytot1 * (ye * df_dt * din - sele)
    dsepdz  = -ytot1 * (ye * df_dt * den  + df_t)


    eele    = ye*free + temp * sele
    deepdt  = temp * dsepdt
    deepdd  = x * df_d + temp * dsepdd
    deepda  = -ye * ytot1 * (free +  df_d * din) + temp * dsepda
    deepdz  = ytot1* (free + ye * df_d * den) + temp * dsepdz




    ! coulomb section:

    ! uniform background corrections only
    ! from yakovlev & shalybkov 1989
    ! lami is the average ion seperation
    ! plasg is the plasma coupling parameter

    z        = forth * pi
    s        = z * xni
    dsdd     = z * dxnidd
    dsda     = z * dxnida

    lami     = 1.0d0/s**third
    inv_lami = 1.0d0/lami
    z        = -third * lami
    lamidd   = z * dsdd/s
    lamida   = z * dsda/s

    plasg    = zbar*zbar*esqu*ktinv*inv_lami
    z        = -plasg * inv_lami
    plasgdd  = z * lamidd
    plasgda  = z * lamida
    plasgdt  = -plasg*ktinv * kerg
    plasgdz  = 2.0d0 * plasg/zbar


    ! ! yakovlev & shalybkov 1989 equations 82, 85, 86, 87
    ! if (plasg .ge. 1.0) then
    !    x        = plasg**(0.25d0)
    !    y        = avo * ytot1 * kerg
    !    ecoul    = y * temp * (a1*plasg + b1*x + c1/x + d1)
    !    pcoul    = third * den * ecoul
    !    scoul    = -y * (3.0d0*b1*x - 5.0d0*c1/x &
    !         + d1 * (log(plasg) - 1.0d0) - e1)

    !    y        = avo*ytot1*kt*(a1 + 0.25d0/plasg*(b1*x - c1/x))
    !    decouldd = y * plasgdd
    !    decouldt = y * plasgdt + ecoul/temp
    !    decoulda = y * plasgda - ecoul/abar
    !    decouldz = y * plasgdz

    !    y        = third * den
    !    dpcouldd = third * ecoul + y*decouldd
    !    dpcouldt = y * decouldt
    !    dpcoulda = y * decoulda
    !    dpcouldz = y * decouldz


    !    y        = -avo*kerg/(abar*plasg)*(0.75d0*b1*x+1.25d0*c1/x+d1)
    !    dscouldd = y * plasgdd
    !    dscouldt = y * plasgdt
    !    dscoulda = y * plasgda - scoul/abar
    !    dscouldz = y * plasgdz


    !    ! yakovlev & shalybkov 1989 equations 102, 103, 104
    ! else if (plasg .lt. 1.0) then
    !    x        = plasg*sqrt(plasg)
    !    y        = plasg**b2
    !    z        = c2 * x - third * a2 * y
    !    pcoul    = -pion * z
    !    ecoul    = 3.0d0 * pcoul/den
    !    scoul    = -avo/abar*kerg*(c2*x -a2*(b2-1.0d0)/b2*y)

    !    s        = 1.5d0*c2*x/plasg - third*a2*b2*y/plasg
    !    dpcouldd = -dpiondd*z - pion*s*plasgdd
    !    dpcouldt = -dpiondt*z - pion*s*plasgdt
    !    dpcoulda = -dpionda*z - pion*s*plasgda
    !    dpcouldz = -dpiondz*z - pion*s*plasgdz

    !    s        = 3.0d0/den
    !    decouldd = s * dpcouldd - ecoul/den
    !    decouldt = s * dpcouldt
    !    decoulda = s * dpcoulda
    !    decouldz = s * dpcouldz

    !    s        = -avo*kerg/(abar*plasg)*(1.5d0*c2*x-a2*(b2-1.0d0)*y)
    !    dscouldd = s * plasgdd
    !    dscouldt = s * plasgdt
    !    dscoulda = s * plasgda - scoul/abar
    !    dscouldz = s * plasgdz
    ! end if


    ! ! bomb proof
    ! x   = prad + pion + pele + pcoul
    ! y   = erad + eion + eele + ecoul
    ! z   = srad + sion + sele + scoul

    ! !        write(6,*) x,y,z
    ! !        if (x .le. 0.0 .or. y .le. 0.0 .or. z .le. 0.0) then
    ! if (x .le. 0.0 .or. y .le. 0.0) then
    !    !        if (x .le. 0.0) then

    !    !         write(6,*)
    !    !         write(6,*) 'coulomb corrections are causing a negative pressure'
    !    !         write(6,*) 'setting all coulomb corrections to zero'
    !    !         write(6,*)

    !    pcoul    = 0.0d0
    !    dpcouldd = 0.0d0
    !    dpcouldt = 0.0d0
    !    dpcoulda = 0.0d0
    !    dpcouldz = 0.0d0
    !    ecoul    = 0.0d0
    !    decouldd = 0.0d0
    !    decouldt = 0.0d0
    !    decoulda = 0.0d0
    !    decouldz = 0.0d0
    !    scoul    = 0.0d0
    !    dscouldd = 0.0d0
    !    dscouldt = 0.0d0
    !    dscoulda = 0.0d0
    !    dscouldz = 0.0d0
    ! end if

    pcoul    = 0.0d0
    dpcouldd = 0.0d0
    dpcouldt = 0.0d0
    dpcoulda = 0.0d0
    dpcouldz = 0.0d0
    ecoul    = 0.0d0
    decouldd = 0.0d0
    decouldt = 0.0d0
    decoulda = 0.0d0
    decouldz = 0.0d0
    scoul    = 0.0d0
    dscouldd = 0.0d0
    dscouldt = 0.0d0
    dscoulda = 0.0d0
    dscouldz = 0.0d0

    ! sum all the gas components
    pgas    = pion + pele + pcoul
    egas    = eion + eele + ecoul
    sgas    = sion + sele + scoul

    dpgasdd = dpiondd + dpepdd + dpcouldd
    dpgasdt = dpiondt + dpepdt + dpcouldt
    dpgasda = dpionda + dpepda + dpcoulda
    dpgasdz = dpiondz + dpepdz + dpcouldz

    degasdd = deiondd + deepdd + decouldd
    degasdt = deiondt + deepdt + decouldt
    degasda = deionda + deepda + decoulda
    degasdz = deiondz + deepdz + decouldz

    dsgasdd = dsiondd + dsepdd + dscouldd
    dsgasdt = dsiondt + dsepdt + dscouldt
    dsgasda = dsionda + dsepda + dscoulda
    dsgasdz = dsiondz + dsepdz + dscouldz




    ! add in radiation to get the total
    pres    = prad + pgas
    ener    = erad + egas
    entr    = srad + sgas

    dpresdd = dpraddd + dpgasdd
    dpresdt = dpraddt + dpgasdt
    dpresda = dpradda + dpgasda
    dpresdz = dpraddz + dpgasdz

    denerdd = deraddd + degasdd
    denerdt = deraddt + degasdt
    denerda = deradda + degasda
    denerdz = deraddz + degasdz

    dentrdd = dsraddd + dsgasdd
    dentrdt = dsraddt + dsgasdt
    dentrda = dsradda + dsgasda
    dentrdz = dsraddz + dsgasdz


    ! for the gas
    ! the temperature and density exponents (c&g 9.81 9.82)
    ! the specific heat at constant volume (c&g 9.92)
    ! the third adiabatic exponent (c&g 9.93)
    ! the first adiabatic exponent (c&g 9.97)
    ! the second adiabatic exponent (c&g 9.105)
    ! the specific heat at constant pressure (c&g 9.98)
    ! and relativistic formula for the sound speed (c&g 14.29)

    zz        = pgas*deni
    zzi       = den/pgas
    chit_gas  = temp/pgas * dpgasdt
    chid_gas  = dpgasdd*zzi
    cv_gas    = degasdt
    x         = zz * chit_gas/(temp * cv_gas)
    gam3_gas  = x + 1.0d0
    gam1_gas  = chit_gas*x + chid_gas
    nabad_gas = x/gam1_gas
    gam2_gas  = 1.0d0/(1.0d0 - nabad_gas)
    cp_gas    = cv_gas * gam1_gas/chid_gas
    z         = 1.0d0 + (egas + light2)*zzi
    sound_gas = clight * sqrt(gam1_gas/z)



    ! for the totals
    zz    = pres*deni
    zzi   = den/pres
    chit  = temp/pres * dpresdt
    chid  = dpresdd*zzi
    cv    = denerdt
    x     = zz * chit/(temp * cv)
    gam3  = x + 1.0d0
    gam1  = chit*x + chid
    nabad = x/gam1
    gam2  = 1.0d0/(1.0d0 - nabad)
    cp    = cv * gam1/chid
    z     = 1.0d0 + (ener + light2)*zzi
    sound = clight * sqrt(gam1/z)



    ! maxwell relations; each is zero if the consistency is perfect
    x   = den * den

    dse = temp*dentrdt/denerdt - 1.0d0

    dpe = (denerdd*x + temp*dpresdt)/pres - 1.0d0

    dsp = -dentrdd*x/dpresdt - 1.0d0
    
    eta_out = etaele
    
    ! ! store this row
    ! ptot_row(j)   = pres
    ! dpt_row(j)    = dpresdt
    ! dpd_row(j)    = dpresdd
    ! dpa_row(j)    = dpresda
    ! dpz_row(j)    = dpresdz
    
    ! etot_row(j)   = ener
    ! det_row(j)    = denerdt
    ! ded_row(j)    = denerdd
    ! dea_row(j)    = denerda
    ! dez_row(j)    = denerdz
    
    ! stot_row(j)   = entr
    ! dst_row(j)    = dentrdt
    ! dsd_row(j)    = dentrdd
    ! dsa_row(j)    = dentrda
    ! dsz_row(j)    = dentrdz
    
    
    ! pgas_row(j)   = pgas
    ! dpgast_row(j) = dpgasdt
    ! dpgasd_row(j) = dpgasdd
    ! dpgasa_row(j) = dpgasda
    ! dpgasz_row(j) = dpgasdz
    
    ! egas_row(j)   = egas
    ! degast_row(j) = degasdt
    ! degasd_row(j) = degasdd
    ! degasa_row(j) = degasda
    ! degasz_row(j) = degasdz
    
    ! sgas_row(j)   = sgas
    ! dsgast_row(j) = dsgasdt
    ! dsgasd_row(j) = dsgasdd
    ! dsgasa_row(j) = dsgasda
    ! dsgasz_row(j) = dsgasdz
    
    
    ! prad_row(j)   = prad
    ! dpradt_row(j) = dpraddt
    ! dpradd_row(j) = dpraddd
    ! dprada_row(j) = dpradda
    ! dpradz_row(j) = dpraddz
    
    ! erad_row(j)   = erad
    ! deradt_row(j) = deraddt
    ! deradd_row(j) = deraddd
    ! derada_row(j) = deradda
    ! deradz_row(j) = deraddz
    
    ! srad_row(j)   = srad
    ! dsradt_row(j) = dsraddt
    ! dsradd_row(j) = dsraddd
    ! dsrada_row(j) = dsradda
    ! dsradz_row(j) = dsraddz
    
    
    ! pion_row(j)   = pion
    ! dpiont_row(j) = dpiondt
    ! dpiond_row(j) = dpiondd
    ! dpiona_row(j) = dpionda
    ! dpionz_row(j) = dpiondz
    
    ! eion_row(j)   = eion
    ! deiont_row(j) = deiondt
    ! deiond_row(j) = deiondd
    ! deiona_row(j) = deionda
    ! deionz_row(j) = deiondz
    
    ! sion_row(j)   = sion
    ! dsiont_row(j) = dsiondt
    ! dsiond_row(j) = dsiondd
    ! dsiona_row(j) = dsionda
    ! dsionz_row(j) = dsiondz
    
    ! xni_row(j)    = xni
    
    ! pele_row(j)   = pele
    ! ppos_row(j)   = 0.0d0
    ! dpept_row(j)  = dpepdt
    ! dpepd_row(j)  = dpepdd
    ! dpepa_row(j)  = dpepda
    ! dpepz_row(j)  = dpepdz
    
    ! eele_row(j)   = eele
    ! epos_row(j)   = 0.0d0
    ! deept_row(j)  = deepdt
    ! deepd_row(j)  = deepdd
    ! deepa_row(j)  = deepda
    ! deepz_row(j)  = deepdz
    
    ! sele_row(j)   = sele
    ! spos_row(j)   = 0.0d0
    ! dsept_row(j)  = dsepdt
    ! dsepd_row(j)  = dsepdd
    ! dsepa_row(j)  = dsepda
    ! dsepz_row(j)  = dsepdz
    
    ! xnem_row(j)   = xnem
    ! xne_row(j)    = xnefer
    ! dxnet_row(j)  = dxnedt
    ! dxned_row(j)  = dxnedd
    ! dxnea_row(j)  = dxneda
    ! dxnez_row(j)  = dxnedz
    ! xnp_row(j)    = 0.0d0
    ! zeff_row(j)   = zbar
    
    ! etaele_row(j) = etaele
    ! detat_row(j)  = detadt
    ! detad_row(j)  = detadd
    ! detaa_row(j)  = detada
    ! detaz_row(j)  = detadz
    ! etapos_row(j) = 0.0d0
    
    ! pcou_row(j)   = pcoul
    ! dpcout_row(j) = dpcouldt
    ! dpcoud_row(j) = dpcouldd
    ! dpcoua_row(j) = dpcoulda
    ! dpcouz_row(j) = dpcouldz
    
    ! ecou_row(j)   = ecoul
    ! decout_row(j) = decouldt
    ! decoud_row(j) = decouldd
    ! decoua_row(j) = decoulda
    ! decouz_row(j) = decouldz
    
    ! scou_row(j)   = scoul
    ! dscout_row(j) = dscouldt
    ! dscoud_row(j) = dscouldd
    ! dscoua_row(j) = dscoulda
    ! dscouz_row(j) = dscouldz
    
    ! plasg_row(j)  = plasg
    
    ! dse_row(j)    = dse
    ! dpe_row(j)    = dpe
    ! dsp_row(j)    = dsp
    
    ! cv_gas_row(j)    = cv_gas
    ! cp_gas_row(j)    = cp_gas
    ! gam1_gas_row(j)  = gam1_gas
    ! gam2_gas_row(j)  = gam2_gas
    ! gam3_gas_row(j)  = gam3_gas
    ! nabad_gas_row(j) = nabad_gas
    ! cs_gas_row(j)    = sound_gas
    
    ! cv_row(j)     = cv
    ! cp_row(j)     = cp
    ! gam1_row(j)   = gam1
    ! gam2_row(j)   = gam2
    ! gam3_row(j)   = gam3
    ! nabad_row(j)  = nabad
    ! cs_row(j)     = sound
    
  end subroutine eos_get_misc

  subroutine eos_get_temp_from_pres(pres, rho, ye, ytot, temp, temp_lo, err_return, return_code_T)
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    ! use module_eos
    ! use mpipara,only: myrank
    ! use unit,only: rho_uni, v_uni, tem_uni
    implicit none
    ! pressure (dyn/cm^2)
    real(8),intent(in) :: pres
    ! density (g/cm^3)
    real(8),intent(in) :: rho
    ! Ye
    real(8),intent(in) :: ye
    ! Ytot = 1/<A>
    real(8),intent(in) :: ytot
    ! T (K)
    real(8),intent(inout) :: temp
    ! Minimum T searched (K)
    real(8),intent(in) :: temp_lo
    ! err
    real(8),intent(out) :: err_return
    ! return code, 0: success, 10: not bracketed, 1: clamped to T_min, 2:likely f=0
    integer,intent(out) :: return_code_T

    integer,parameter :: itrlim=20
    real(8),parameter :: tol=1d-10
    
    real(8) :: logtemp
    
    
    real(8) :: pres_min, pres_lo
    real(8) :: eps,cs2,mexc,entr
!!!

    ! As only pressure is used, mexc is just a dummy argument
    mexc = 0d0
    
    call eos_all(rho, temp_lo, ye, ytot, mexc, &
         eps, pres_lo, cs2, entr)

    if(pres < pres_lo)then
       ! write(6,*) "eps < eps_min", eps, eps_min
       ! temp = 10d0**logtemp_lo
       temp = temp_lo
       err_return = 0d0
       return_code_T = 1 ! clamped to T_min
       return
    endif

    call eos_all(rho, temp_eos_min, ye, ytot, mexc, &
         eps, pres_min, cs2, entr)

    if (pres_lo <= pres_min) then
       write(6,*) 'ERROR: invalid pres->T inversion floor'
       write(6,*) 'rho, temp_lo, pres_lo, pres_min = ', &
            rho, temp_lo, pres_lo, pres_min
       error stop
    endif

    block
      
      integer :: itr
      real(8) :: x, err, x_min, x_max, itr_out
      integer :: side
      real(8) :: x1,x2, f1,f2,f,f_max, tol_bar

      x_min = log10(temp_lo)
      x_max = log10(temp_eos_max)

      side = 0
      x1=x_min
      x2=x_max
      
      logtemp = x1
      call master_function_temp_from_pres(logtemp, f, pres, rho, ye, ytot, mexc, pres_min)
      f1 = f

      logtemp = x2
      call master_function_temp_from_pres(logtemp, f, pres, rho, ye, ytot, mexc, pres_min)
      f2 = f
      
      f_max = max(abs(f1),abs(f2))
      
      if(f1*f2>0d0)then
         
         if(abs(f1) < abs(f2))then
            f=f1
            x=x1
         else
            f=f2
            x=x2
         endif
         write(6,*) "not bracketed (pres->T)", f1,f2
         itr_out = 0
         tol_bar = tol
         do itr=1,10
            tol_bar = tol_bar*10d0
            !write(6,*) itr,tol_bar,abs(f)
            if(abs(f) <= tol_bar)then
               err = tol_bar
            endif
         enddo
         temp = 10.d0**x
         err_return = abs(f)
         return_code_T = 10         
         return
      else
         
         loop_itr: do itr=1,itrlim
              
            x = (x1*f2 - x2*f1)/(f2-f1)
            
            logtemp = x
            call master_function_temp_from_pres(logtemp, f, pres, rho, ye, ytot, mexc, pres_min)
            
            ! write(6,'(i5,99es12.4)') itr, f, f1, f2
            
            if (.not. ieee_is_finite(f)) then
               write(*,*) 'ERROR: non-finite f in pres->T inversion'
               write(*,*) 'x, x1, x2 = ', x, x1, x2
               write(*,*) 'f, f1, f2 = ', f, f1, f2
               temp = temp_lo
               err_return = huge(1d0)
               return_code_T = 20
               return
            endif

            if(f*f1>0d0)then
               x1 = x
               f1 = f
               if(side==1) f2 = 0.5d0*f2
               side = 1
            elseif(f*f2>0d0)then
               x2 = x
               f2 = f
               if(side==2) f1 = 0.5d0*f1
               side = 2
            else
               ! write(6,*) "already root?"
               ! write(6,*) f,f1,f2
               temp = 10d0**logtemp
               err_return = 0d0
               return_code_T = 0
               return
            endif

            ! write(6,'(i5,99es12.4)') itr, x, abs(x1-x2)/abs(x1+x2)
            
            if( abs(x1-x2) <= tol*abs(x1+x2) ) then
               err = abs(x1-x2)/abs(x1+x2)
               return_code_T = 0

               logtemp = x
               temp = 10d0**logtemp
               err_return = err
               return
            endif
            itr_out = itr
         enddo loop_itr
         if(itr>itrlim) err = abs(x1-x2)/abs(x1+x2)
         ! write(6,*) err
         err = abs(x1-x2)/abs(x1+x2)
         logtemp = x
         temp = 10d0**logtemp
         err_return = err
         return_code_T = 2 ! Reached iteration limit
         return
         
      endif

    end block
    
  end subroutine eos_get_temp_from_pres

  subroutine master_function_temp_from_pres(logtemp, f, pres_target, rho, ye, ytot, mexc, pres_min)
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    ! log10(T (K))
    real(8),intent(in) :: logtemp
    ! master function
    real(8),intent(out) :: f
    ! target pressure (dyn/cm^2)
    real(8),intent(in) :: pres_target
    ! rho (g/cm^3)
    real(8),intent(in) :: rho
    ! Ye
    real(8),intent(in) :: ye
    ! Ytot = 1/<A>
    real(8),intent(in) :: ytot
    ! Mass excess per baryon (MeV/b)
    real(8),intent(in) :: mexc
    ! minimum pressure (dyn/cm^2)
    real(8),intent(in) :: pres_min

    real(8) :: temp_in, rho_in, eps_dummy, pres, cs2_dummy, entr_dummy
    
    temp_in = 10d0**logtemp
    rho_in = rho
    call eos_all(rho_in, temp_in, ye, ytot, mexc, &
         eps_dummy, pres, cs2_dummy, entr_dummy)
    
    block    
      real(8) :: dp_target, dp
      real(8) :: pscale, p_floor

      pscale  = max(abs(pres_target), abs(pres_min), 1d0)
      p_floor = 32d0 * spacing(pscale)

      dp_target = pres_target - pres_min
      dp        = pres        - pres_min

      if (.not. ieee_is_finite(dp)) then
         write(6,*) 'ERROR: non-finite pressure difference'
         error stop
      endif

      ! 数 ulp 程度の負値は丸め誤差として扱う
      if (dp < -p_floor) then
         write(6,*) 'ERROR: pressure substantially below pres_min'
         write(6,'(a,es26.17e3)') 'dp      = ', dp
         write(6,'(a,es26.17e3)') 'p_floor = ', p_floor
         error stop
      endif

      dp_target = max(dp_target, 0d0)
      dp        = max(dp,        0d0)

      f = log10(dp_target + p_floor) - log10(dp + p_floor)
    end block
    ! f = log10( max(pres_target,pres_min)-pres_min+1d-10)-log10(pres-pres_min+1d-10)
    
  end subroutine master_function_temp_from_pres


  subroutine eos_get_temp_from_eps(eps, rho, ye, ytot, mexc, temp, temp_lo, err_return, return_code_T)
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    ! epsilon (erg/g)
    real(8),intent(in) :: eps
    ! density [g cm^-3]
    real(8),intent(in) :: rho
    ! Ye
    real(8),intent(in) :: ye
    ! Ytot = 1/<A>
    real(8),intent(in) :: ytot
    ! mexc (MeV)
    real(8),intent(in) :: mexc
    ! T (K)
    real(8),intent(inout) :: temp
    ! Minimum T searched (K)
    real(8),intent(in) :: temp_lo
    ! err
    real(8),intent(out) :: err_return
    ! return code, 0: success, 10: not bracketed, 1: clamped to T_min, 2:likely f=0
    integer,intent(out) :: return_code_T
    
    integer,parameter :: itrlim=20
    real(8),parameter :: tol=1d-10
    ! integer :: ie,ied0
    real(8) :: logtemp,logtemp_guess
    real(8) :: rho_cgs

    real(8) :: eps_min, eps_lo
    
    rho_cgs = rho

    block
      real(8) :: pres,cs2,entr
      call eos_all(rho_cgs, temp_lo, ye, ytot, mexc, &
           eps_lo, pres, cs2, entr)
    end block

    if(eps < eps_lo)then
       ! write(6,*) "eps < eps_lo", eps, eps_lo
       ! temp = 10d0**logtemp_lo
       temp = temp_lo
       err_return = 0d0
       return_code_T = 1 ! clamped to T_min
       return
    endif

    block
      real(8) :: pres,cs2,entr
      call eos_all(rho_cgs, temp_eos_min, ye, ytot, mexc, &
           eps_min, pres, cs2, entr)
    end block

    if (eps_lo <= eps_min) then
       write(6,*) 'ERROR: invalid eps->T inversion floor'
       write(6,*) 'rho, temp_lo, eps_lo, eps_min = ', &
            rho, temp_lo, eps_lo, eps_min
       error stop
    endif

    block
      
      integer :: itr
      real(8) :: x, err, x_min, x_max, itr_out
      integer :: side
      real(8) :: x1,x2, f1,f2,f,f_max, tol_bar

      x_min = log10(temp_lo)
      x_max = log10(temp_eos_max)

      side = 0
      x1=x_min
      x2=x_max
      
      logtemp = x1
      call master_function_temp_from_eps(logtemp, f, eps, rho, ye, ytot, mexc, eps_min)
      f1 = f

      logtemp = x2
      call master_function_temp_from_eps(logtemp, f, eps, rho, ye, ytot, mexc, eps_min)
      f2 = f

      ! write(6,'("T_min,T_max,f1,f2=",99es12.4)') x1,x2,f1,f2
      
      f_max = max(abs(f1),abs(f2))
     
      if(f1*f2>0d0)then
         
         if(abs(f1) < abs(f2))then
            f=f1
            x=x1
         else
            f=f2
            x=x2
         endif
         write(6,*) "not bracketed (eps->T)", f1,f2
         itr_out = 0
         tol_bar = tol
         do itr=1,10
            tol_bar = tol_bar*10d0
            !write(6,*) itr,tol_bar,abs(f)
            if(abs(f) <= tol_bar)then
               err = tol_bar
            endif
         enddo
         temp = 10.d0**x
         err_return = abs(f)
         return_code_T = 10         
         return
      else
         
         loop_itr: do itr=1,itrlim
              
            x = (x1*f2 - x2*f1)/(f2-f1)
            
            logtemp = x
            call master_function_temp_from_eps(logtemp, f, eps, rho, ye, ytot, mexc, eps_min)
            
            if (.not. ieee_is_finite(f)) then
               write(*,*) 'ERROR: non-finite f in eps->T inversion'
               write(*,*) 'x, x1, x2 = ', x, x1, x2
               write(*,*) 'f, f1, f2 = ', f, f1, f2
               temp = temp_lo
               err_return = huge(1d0)
               return_code_T = 20
               return
            endif
            ! write(6,'(i5,99es12.4)') itr, f, f1, f2
            
            if(f*f1>0d0)then
               x1 = x
               f1 = f
               if(side==1) f2 = 0.5d0*f2
               side = 1
            elseif(f*f2>0d0)then
               x2 = x
               f2 = f
               if(side==2) f1 = 0.5d0*f1
               side = 2
            else
               ! write(6,*) "already root?"
               ! write(6,*) f,f1,f2
               temp = 10d0**logtemp
               err_return = 0d0
               return_code_T = 0
               return
            endif

            ! write(6,'(i5,99es12.4)') itr, x, abs(x1-x2)/abs(x1+x2)
            
            if( abs(x1-x2) <= tol*abs(x1+x2) ) then
               err = abs(x1-x2)/abs(x1+x2)
               return_code_T = 0

               logtemp = x
               temp = 10d0**logtemp
               err_return = err
               return
            endif
            itr_out = itr
         enddo loop_itr
        
         err = abs(x1-x2)/abs(x1+x2)
         logtemp = x
         temp = 10d0**logtemp
         err_return = err
         return_code_T = 2 ! Reached iteration limit
         return
      endif
    
    end block
    
    ! block
    !   integer,parameter :: n = 10
    !   integer :: i
    !   do i=1,n
    !      logtemp = logtemp_lo + (logtemp_max-logtemp_lo)*dble(i-1)/dble(n-1)
    !      write(6,*) logtemp, func(logtemp)
    !   enddo
    ! end block
    ! stop

    ! call illinois_rootfinding(logtemp,err,itr_tem,logtemp_lo,logtemp_max,tol,itrlim_tem,return_code,func)
    


    ! write(6,'("eps,rho -> temp:",i4,99es12.4)') itr_tem,logtemp,temp,err
    
    
  end subroutine eos_get_temp_from_eps

  subroutine master_function_temp_from_eps(logtemp, f, eps_target, rho, ye, ytot, mexc, eps_min)
    
    ! log10(T (K))
    real(8),intent(in) :: logtemp
    ! master function
    real(8),intent(out) :: f
    ! target eps (c^2)
    real(8),intent(in) :: eps_target
    ! rho (g cm^-3)
    real(8),intent(in) :: rho
    ! Ye
    real(8),intent(in) :: ye
    ! Ytot = 1/<A>
    real(8),intent(in) :: ytot
    ! mean mass excess (MeV)
    real(8),intent(in) :: mexc
    ! minimum epsilon (c^2)
    real(8),intent(in) :: eps_min
    
    
    real(8) :: temp_k, rho_cgs, eps, pres_dummy, cs2_dummy, entr_dummy
    
    
    temp_k = 10d0**logtemp
    rho_cgs = rho
    call eos_all(rho_cgs, temp_k, ye, ytot, mexc, &
         eps, pres_dummy, cs2_dummy, entr_dummy)
    
    block
      use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    
      real(8) :: escale, e_floor, de_target, de
      escale  = max(abs(eps_target), abs(eps_min), 1d0)
      e_floor = 32d0 * spacing(escale)
      
      de_target = eps_target - eps_min
      de        = eps        - eps_min
      
      if (.not. ieee_is_finite(de)) then
         write(6,*) 'ERROR: non-finite energy difference'
         error stop
      endif
      
      ! A small negative value is allowed as roundoff.
      if (de < -e_floor) then
         write(6,*) 'ERROR: energy substantially below eps_min'
         write(6,'(a,es26.17e3)') 'logtemp    = ', logtemp
         write(6,'(a,es26.17e3)') 'temp       = ', 10.d0**logtemp
         write(6,'(a,es26.17e3)') 'rho        = ', rho
         write(6,'(a,es26.17e3)') 'ye         = ', ye
         write(6,'(a,es26.17e3)') 'ytot       = ', ytot
         write(6,'(a,es26.17e3)') 'mexc       = ', mexc
         write(6,'(a,es26.17e3)') 'eps_target = ', eps_target
         write(6,'(a,es26.17e3)') 'eps        = ', eps
         write(6,'(a,es26.17e3)') 'eps_min    = ', eps_min
         write(6,'(a,es26.17e3)') 'de         = ', de
         write(6,'(a,es26.17e3)') 'e_floor    = ', e_floor
         error stop
      endif
      
      de_target = max(de_target, 0d0)
      de        = max(de,        0d0)
      
      f = log10(de_target + e_floor) - log10(de        + e_floor)

    end block
    !f = log10(eps) - log10(eps_target)
    ! write(6,'("T, rho, eps_target, eps, eps_min = ",99es12.4)') logtemp, rho_cgs, ye, ytot, mexc, eps_target, eps, eps_min
    
    ! f = log10( max(eps_target,eps_min)-eps_min+1d-10)-log10(eps-eps_min+1d-10)

  end subroutine master_function_temp_from_eps


  function eos_get_eps_at_temp_min(rho,ye,ytot,mexc) result(eps)
    ! rho [g cm^-3]
    real(8) :: rho
    ! Ye, Ytot, mass excess per baryon [MeV]
    real(8),intent(in) :: ye, ytot, mexc
    
    ! specific internal energy [erg g^-1]
    real(8) :: eps

    real(8) :: eps_out, pres_out, cs2_out, entr_out
        
    call eos_all(rho, temp_eos_min, ye, ytot, mexc, &
         eps_out, pres_out, cs2_out, entr_out)
    
    eps = eps_out
    
  end function eos_get_eps_at_temp_min

end module module_eos_helmholtz
