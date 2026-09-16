module module_stat_weight_HS

  implicit none
  private

  public :: get_stat_weight_HS

contains

  subroutine get_stat_weight_HS(t9, a, z, n, mexc, g)

    use const, only: mev2k, mpmev, mnmev, mumev, pi

    real(8), intent(in)  :: t9
    real(8), intent(in)  :: a, z, n, mexc
    real(8), intent(out) :: g

    real(8), parameter :: c1 = 0.2d0
    real(8), parameter :: c2 = 0.8d0

    real(8) :: temp_mev
    real(8) :: aa, emax, bind, g0, iexc
    integer :: ia, iz

    ia = nint(a)
    iz = nint(z)

    temp_mev = t9*1d9/mev2k

    ! free neutron / proton
    if (ia == 1) then
       g = 2d0
       return
    endif

    ! Fai-Randrup / HS ground-state prescription
    if (mod(ia,2) == 0) then
       g0 = 1d0
    else
       g0 = 2d0
    endif

    ! deuteron
    if (ia == 2 .and. iz == 1) then
       g0 = 3d0
    endif

    ! nuclear binding energy [MeV]
    bind = z*mpmev + n*mnmev &
         - (a*mumev + mexc)

    emax = max(bind, 0d0)

    if (temp_mev <= 0d0 .or. emax <= 0d0) then
       g = g0
       return
    endif

    aa = a/8d0 * (1d0 - c2*a**(-1d0/3d0))

    iexc = excited_HS(temp_mev, aa, emax)

    g = g0 + c1/a**(5d0/3d0)*iexc

  end subroutine get_stat_weight_HS


  function excited_HS(temp, aa, emax) result(val)

    use const, only: pi

    real(8), intent(in) :: temp, aa, emax
    real(8) :: val

    real(8) :: cc, s0, s1
    real(8) :: derf, term1, term2, xexp

    cc = temp*sqrt(aa/2d0)

    s0 = -sqrt(aa*temp/2d0)
    s1 = (sqrt(emax)-cc)/sqrt(temp)

    if (s0 < 0d0 .and. s1 < 0d0) then
       derf = erfc(-s1) - erfc(-s0)
    elseif (s0 > 0d0 .and. s1 > 0d0) then
       derf = erfc(s0) - erfc(s1)
    else
       derf = erf(s1) - erf(s0)
    endif
    
    term1 = cc*sqrt(pi*temp)*derf
    term2 = temp*(exp(-s0*s0)-exp(-s1*s1))

    xexp = aa*temp/2.d0
    
    if (xexp > log(huge(1.d0)) - 10.d0) then
       write(*,*) "ERROR: HS partition function exponent too large", &
            aa, temp, xexp
       error stop
    endif

    val = exp(aa*temp/2.d0) * (term1 + term2)

    if(val<=0d0)then
       write(6,*) "val<0.0", temp, aa, emax, s0, s1, term1, term2, term1+term2
       error stop
    endif

  end function excited_HS

end module module_stat_weight_HS
