program make_compose_helmholtz
  use module_nse, only: nse_network_t, calc_nse, nse_init_winvne, calc_nse_nested_1d
  use module_nuclear_data_winvne, only: init_winvne
  use module_ptf_rauscher, only: init_ptf_rauscher
  use module_nuclear_data_HS, only: init_nuclear_data_HS
  use module_eos_helmholtz, only: init_eos, eos_all, eos_get_misc
  use module_stat_weight_policy, only: stat_weight_policy_t, &
       STAT_WEIGHT_HS, STAT_WEIGHT_WINVNE
  use module_nuclear_mass_policy, only: nuclear_mass_policy_t, &
       NUCLEAR_MASS_NONE, NUCLEAR_MASS_WINVNE
  use module_nse_species_policy, only: nse_species_policy_t, &
       NSE_SPECIES_LEGACY
  use module_compose_hdf5, only: write_compose_hdf5
  use const, only: mu, mev2erg, mev2k, clight, mumev, mnmev, mpmev, memev
  implicit none

  integer, parameter :: itrlim = 300
  real(8), parameter :: tol = 1.d-10
  real(8), parameter :: norm_tol = 1.d-10

  character(512) :: fn_para
  character(512) :: fn_winv, fn_rauscher, fn_hs, fn_helm, fn_out

  integer :: nnb, nt, nyq
  real(8) :: nb_min, nb_max, t_min, t_max, yq_min, yq_max

  real(8), allocatable :: nb(:), t(:), yq(:)
  real(8), allocatable :: q1(:,:,:), q2(:,:,:), q3(:,:,:), q4(:,:,:)
  real(8), allocatable :: q5(:,:,:), q6(:,:,:), q7(:,:,:), cs2(:,:,:)
  real(8), allocatable :: ye_tab(:,:,:), yn(:,:,:), yp(:,:,:)
  real(8), allocatable :: yh2(:,:,:), yh3(:,:,:), yhe3(:,:,:), yhe4(:,:,:)
  real(8), allocatable :: ynuc(:,:,:), anuc(:,:,:), znuc(:,:,:), abar(:,:,:)
  real(8), allocatable :: mue(:,:,:)

  real(8), allocatable :: xnse(:)
  type(nse_network_t) :: net
  type(stat_weight_policy_t) :: stat_weight_policy
  type(nuclear_mass_policy_t) :: nuclear_mass_policy
  type(nse_species_policy_t) :: species_policy

  integer :: inb, it, iyq
  integer :: itr_out
  integer :: fallback_count
  logical :: nsefail
  real(8) :: err_out, xn_guess, xp_guess, xn_out, xp_out
  real(8) :: rho, temp_k, ye
  real(8) :: mexc, ytot
  real(8) :: eps, pres, cs2_cgs, entr, eta_e
  real(8) :: baryon_sum, charge_sum
  real(8) :: max_baryon_err, max_charge_err

  integer :: u_summary

  call get_command_argument(1, fn_para)
  if (len_trim(fn_para) == 0) then
     write(*,*) "usage: make_compose_helmholtz <parameter-file>"
     error stop
  endif

  call read_parameter_file(trim(fn_para), &
       fn_winv, fn_rauscher, fn_hs, fn_helm, fn_out, &
       nnb, nb_min, nb_max, nt, t_min, t_max, nyq, yq_min, yq_max)

  call validate_grid_parameters(nnb, nb_min, nb_max, &
       nt, t_min, t_max, nyq, yq_min, yq_max)

  allocate(nb(nnb), t(nt), yq(nyq))
  call make_log_grid(nb, nb_min, nb_max)
  call make_log_grid(t,  t_min,  t_max)
  call make_linear_grid(yq, yq_min, yq_max)

  ! Match the choices currently used by compose_extraction.f90:
  !   nuclear masses: WinVNE
  !   statistical weights: HS, with WinVNE fallback
  !   species set: legacy large NSE set
  nuclear_mass_policy%primary  = NUCLEAR_MASS_WINVNE
  nuclear_mass_policy%fallback = NUCLEAR_MASS_NONE

  stat_weight_policy%primary  = STAT_WEIGHT_HS
  stat_weight_policy%fallback = STAT_WEIGHT_WINVNE

  species_policy%mode = NSE_SPECIES_LEGACY

  call init_winvne(trim(fn_winv))
  call init_nuclear_data_HS(trim(fn_hs))
  call init_ptf_rauscher(trim(fn_rauscher))
  call nse_init_winvne(net, stat_weight_policy, species_policy, nuclear_mass_policy)
  call init_eos(trim(fn_helm))

  allocate(xnse(net%n_spec))

  allocate(q1(nt,nyq,nnb), q2(nt,nyq,nnb), q3(nt,nyq,nnb), q4(nt,nyq,nnb))
  allocate(q5(nt,nyq,nnb), q6(nt,nyq,nnb), q7(nt,nyq,nnb), cs2(nt,nyq,nnb))
  allocate(ye_tab(nt,nyq,nnb), yn(nt,nyq,nnb), yp(nt,nyq,nnb))
  allocate(yh2(nt,nyq,nnb), yh3(nt,nyq,nnb), yhe3(nt,nyq,nnb), yhe4(nt,nyq,nnb))
  allocate(ynuc(nt,nyq,nnb), anuc(nt,nyq,nnb), znuc(nt,nyq,nnb), abar(nt,nyq,nnb))
  allocate(mue(nt,nyq,nnb))

  q1 = 0.d0; q2 = 0.d0; q3 = 0.d0; q4 = 0.d0
  q5 = 0.d0; q6 = 0.d0; q7 = 0.d0; cs2 = 0.d0
  ye_tab = 0.d0; yn = 0.d0; yp = 0.d0
  yh2 = 0.d0; yh3 = 0.d0; yhe3 = 0.d0; yhe4 = 0.d0
  ynuc = 0.d0; anuc = 0.d0; znuc = 0.d0; abar = 0.d0
  mue = 0.d0

  fallback_count = 0
  max_baryon_err = 0.d0
  max_charge_err = 0.d0

  write(*,'(a,3(i0,1x))') "grid (nb, yq, t) = ", nnb, nyq, nt
  write(*,'(a,2es14.6)') "nb [fm^-3] = ", nb(1), nb(nnb)
  write(*,'(a,2es14.6)') "T  [MeV]   = ", t(1),  t(nt)
  write(*,'(a,2es14.6)') "Yq         = ", yq(1), yq(nyq)
  write(*,'(a,i0)') "NSE species = ", net%n_spec

  open(newunit=u_summary, file= "summary.dat", status="replace", action="write")

  !$omp parallel default(none) &
  !$omp shared(net, u_summary, nyq, nt, nnb, yq, t, nb, ye_tab, q1, q2, q7, q6, cs2, mue, &
  !$omp   yn, yp, yh2, yh3, yhe3, yhe4, ynuc, anuc, znuc, abar) &
  !$omp private(ye, temp_k, rho, xnse, nsefail, itr_out, err_out, xn_out, xp_out, xn_guess, xp_guess, &
  !$omp   mexc, ytot, charge_sum, baryon_sum, max_baryon_err, max_charge_err, eps, pres, cs2_cgs, entr, eta_e)
  !$omp do collapse(2) schedule(dynamic,1)
  do iyq = 1, nyq
     do it = 1, nt
        ye = yq(iyq)
        temp_k = t(it) * mev2k

        do inb = 1, nnb
           rho = mu * nb(inb) * 1.d39

           ! Continue the neutron/proton chemical-potential variables along
           ! the density direction whenever possible.

           if (inb == 1) then

              call calc_nse_nested_1d( &
                   net, rho, temp_k, ye, itrlim, tol, &
                   xnse, nsefail, &
                   itr_out=itr_out, err_out=err_out, &
                   xn_out=xn_out, xp_out=xp_out)
           else
              call calc_nse_nested_1d( &
                   net, rho, temp_k, ye, itrlim, tol, &
                   xnse, nsefail, &
                   itr_out=itr_out, err_out=err_out, &
                   xn_guess=xn_guess, xp_guess=xp_guess, &
                   xn_out=xn_out, xp_out=xp_out)

           endif

           ! if (inb == 1) then
           !    call calc_nse(net, rho, temp_k, ye, itrlim, tol, &
           !         xnse, nsefail, .false., &
           !         itr_out=itr_out, err_out=err_out, &
           !         xn_out=xn_out, xp_out=xp_out)
           ! else
           !    call calc_nse(net, rho, temp_k, ye, itrlim, tol, &
           !         xnse, nsefail, .false., &
           !         itr_out=itr_out, err_out=err_out, &
           !         xn_guess=xn_guess, xp_guess=xp_guess, &
           !         xn_out=xn_out, xp_out=xp_out)
           ! endif

           ! if (nsefail) then
           !    fallback_count = fallback_count + 1
           !    call calc_nse(net, rho, temp_k, ye, itrlim, tol, &
           !         xnse, nsefail, .true., &
           !         itr_out=itr_out, err_out=err_out, &
           !         xn_out=xn_out, xp_out=xp_out)
           ! endif

           if (nsefail) then
              write(*,'(a,3i7,3es18.9)') &
                   "NSE failed at (inb,iyq,it), nb,T,Yq = ", &
                   inb, iyq, it, nb(inb), t(it), ye
              write(*,'(a,es18.9)') "final residual = ", err_out
              error stop
           endif

           xn_guess = xn_out
           xp_guess = xp_out

           call compose_moments(net, xnse, ye, &
                mexc, ytot, &
                yn(it,iyq,inb), yp(it,iyq,inb), &
                yh2(it,iyq,inb), yh3(it,iyq,inb), &
                yhe3(it,iyq,inb), yhe4(it,iyq,inb), &
                ynuc(it,iyq,inb), anuc(it,iyq,inb), &
                znuc(it,iyq,inb), abar(it,iyq,inb), &
                baryon_sum, charge_sum)

           max_baryon_err = max(max_baryon_err, abs(baryon_sum - 1.d0))
           max_charge_err = max(max_charge_err, abs(charge_sum - ye))

           if (abs(baryon_sum - 1.d0) > norm_tol .or. &
                abs(charge_sum - ye) > norm_tol) then
              write(*,'(a,3i7)') "composition check failed at ", inb, iyq, it
              write(*,'(a,2es18.9)') "sum(A_i Y_i), target = ", baryon_sum, 1.d0
              write(*,'(a,2es18.9)') "sum(Z_i Y_i), target = ", charge_sum, ye
              error stop
           endif

           ye_tab(it,iyq,inb) = ye

           call eos_all(rho, temp_k, ye, ytot, mexc, &
                eps, pres, cs2_cgs, entr)
           call eos_get_misc(rho, temp_k, ye, ytot, mexc, eta_e)

           ! PyCompOSE conventions used by the DD2.h5 table.
           ! Q1 = P/n_b [MeV]
           q1(it,iyq,inb) = pres / (nb(inb) * 1.d39) / mev2erg

           ! Q2 = entropy per baryon [k_B]
           q2(it,iyq,inb) = entr

           ! Q7 = E_b / m_n - 1, where eps excludes the m_u c^2 baseline.
           q7(it,iyq,inb) = &
                (mumev + eps * mu / mev2erg) / mnmev - 1.d0

           ! Q6 = F_b / m_n - 1 = Q7 - T S / m_n.
           q6(it,iyq,inb) = q7(it,iyq,inb) &
                - t(it) * q2(it,iyq,inb) / mnmev

           ! eos_all returns c_s^2 in cgs units.
           cs2(it,iyq,inb) = cs2_cgs / clight**2

           ! Timmes Helmholtz eta_e is the kinetic electron chemical
           ! potential divided by kT; add m_e c^2 for the relativistic
           ! chemical potential used by CompOSE.
           mue(it,iyq,inb) = memev + eta_e * t(it)

           write(u_summary,'(99es20.10e3)') nb(inb), t(it), ye, xn_out, xp_out, yn(it,iyq,inb), yp(it,iyq,inb), ynuc(it,iyq,inb)*anuc(it,iyq,inb), anuc(it,iyq,inb), &
                znuc(it,iyq,inb), abar(it,iyq,inb)

        enddo

        write(*,'(a,2i6,a,es12.4,a,es12.4)') &
             "completed iyq,it = ", iyq, it, &
             "  Yq=", ye, "  T[MeV]=", t(it)
     enddo
  enddo
  !$omp end do
  !$omp end parallel

  ! Q5 is obtained from dF_b/dYq at fixed (nb,T), including the
  ! response of the NSE composition.  For charge-neutral matter this is
  ! the electron-lepton chemical potential mu_l.  Then
  !   mu_e = mu_l - mu_q,
  !   F_b  = -P/n_b + mu_b + Yq*mu_l,
  ! which determines Q4 and Q3.
  call fill_chemical_potentials(t, yq, q1, q6, mue, q3, q4, q5)

  call write_compose_hdf5(trim(fn_out), nb, t, yq, mnmev, mpmev, &
       q1, q2, q3, q4, q5, q6, q7, cs2, &
       ye_tab, yn, yp, yh2, yh3, yhe3, yhe4, ynuc, &
       anuc, znuc, abar)

  write(*,'(a,1x,a)') "wrote", trim(fn_out)
  write(*,'(a,i0)') "NSE TNA fallbacks = ", fallback_count
  write(*,'(a,es12.4)') "max |sum A_i Y_i - 1| = ", max_baryon_err
  write(*,'(a,es12.4)') "max |sum Z_i Y_i - Yq| = ", max_charge_err

contains

  subroutine read_parameter_file(fn, &
       fn_winv, fn_rauscher, fn_hs, fn_helm, fn_out, &
       nnb, nb_min, nb_max, nt, t_min, t_max, nyq, yq_min, yq_max)

    character(*), intent(in) :: fn
    character(*), intent(out) :: fn_winv, fn_rauscher, fn_hs, fn_helm, fn_out
    integer, intent(out) :: nnb, nt, nyq
    real(8), intent(out) :: nb_min, nb_max, t_min, t_max, yq_min, yq_max
    real(8) :: lognb_min, lognb_max, logt_min, logt_max

    integer :: iu, ios

    open(newunit=iu, file=fn, status="old", action="read", iostat=ios)
    if (ios /= 0) then
       write(*,*) "cannot open parameter file: ", trim(fn)
       error stop
    endif

    read(iu,*); read(iu,'(a)') fn_winv
    read(iu,*); read(iu,'(a)') fn_rauscher
    read(iu,*); read(iu,'(a)') fn_hs
    read(iu,*); read(iu,'(a)') fn_helm
    read(iu,*); read(iu,'(a)') fn_out
    read(iu,*); read(iu,*) nnb, lognb_min, lognb_max
    read(iu,*); read(iu,*) nt,  logt_min,  logt_max
    read(iu,*); read(iu,*) nyq, yq_min, yq_max

    close(iu)

    nb_min = 10d0**lognb_min
    nb_max = 10d0**lognb_max

    t_min = 10d0**logt_min
    t_max = 10d0**logt_max

    fn_winv     = trim(adjustl(fn_winv))
    fn_rauscher = trim(adjustl(fn_rauscher))
    fn_hs       = trim(adjustl(fn_hs))
    fn_helm     = trim(adjustl(fn_helm))
    fn_out      = trim(adjustl(fn_out))

  end subroutine read_parameter_file


  subroutine validate_grid_parameters(nnb, nb_min, nb_max, &
       nt, t_min, t_max, nyq, yq_min, yq_max)
    integer, intent(in) :: nnb, nt, nyq
    real(8), intent(in) :: nb_min, nb_max, t_min, t_max, yq_min, yq_max

    if (nnb < 1 .or. nt < 1) error stop "nnb and nt must be >= 1"
    if (nyq < 3) error stop "nyq must be >= 3 to construct Q5 with a 3-point derivative"
    if (nb_min <= 0.d0 .or. nb_max < nb_min) error stop "invalid nb range"
    if (t_min <= 0.d0 .or. t_max < t_min) error stop "invalid temperature range"
    if (yq_min <= 0.d0 .or. yq_max >= 1.d0 .or. yq_max <= yq_min) &
         error stop "require 0 < yq_min < yq_max < 1"
  end subroutine validate_grid_parameters


  subroutine make_log_grid(x, xmin, xmax)
    real(8), intent(out) :: x(:)
    real(8), intent(in) :: xmin, xmax
    integer :: i, n

    n = size(x)
    if (n == 1) then
       x(1) = xmin
       return
    endif

    do i = 1, n
       x(i) = 10.d0**(log10(xmin) + &
            (log10(xmax) - log10(xmin))*dble(i-1)/dble(n-1))
    enddo
  end subroutine make_log_grid


  subroutine make_linear_grid(x, xmin, xmax)
    real(8), intent(out) :: x(:)
    real(8), intent(in) :: xmin, xmax
    integer :: i, n

    n = size(x)
    if (n == 1) then
       x(1) = xmin
       return
    endif

    do i = 1, n
       x(i) = xmin + (xmax-xmin)*dble(i-1)/dble(n-1)
    enddo
  end subroutine make_linear_grid


  subroutine compose_moments(net, x, ye_target, &
       mexc, ytot, yn, yp, yh2, yh3, yhe3, yhe4, &
       ynuc, anuc, znuc, abar, baryon_sum, charge_sum)

    type(nse_network_t), intent(in) :: net
    real(8), intent(in) :: x(net%n_spec)
    real(8), intent(in) :: ye_target
    real(8), intent(out) :: mexc, ytot
    real(8), intent(out) :: yn, yp, yh2, yh3, yhe3, yhe4
    real(8), intent(out) :: ynuc, anuc, znuc, abar
    real(8), intent(out) :: baryon_sum, charge_sum

    integer :: k, ia, iz
    real(8) :: yk
    real(8) :: a_num, z_num

    yn = 0.d0; yp = 0.d0
    yh2 = 0.d0; yh3 = 0.d0; yhe3 = 0.d0; yhe4 = 0.d0
    ynuc = 0.d0
    a_num = 0.d0
    z_num = 0.d0

    mexc = 0.d0
    ytot = 0.d0
    baryon_sum = 0.d0
    charge_sum = 0.d0

    do k = 1, net%n_spec
       yk = x(k) / net%a(k)
       ia = nint(net%a(k))
       iz = nint(net%z(k))

       mexc = mexc + net%mexc(k) * yk
       ytot = ytot + yk
       baryon_sum = baryon_sum + net%a(k) * yk
       charge_sum = charge_sum + net%z(k) * yk

       if (ia == 1 .and. iz == 0) then
          yn = yn + yk
       elseif (ia == 1 .and. iz == 1) then
          yp = yp + yk
       elseif (ia == 2 .and. iz == 1) then
          yh2 = yh2 + yk
       elseif (ia == 3 .and. iz == 1) then
          yh3 = yh3 + yk
       elseif (ia == 3 .and. iz == 2) then
          yhe3 = yhe3 + yk
       elseif (ia == 4 .and. iz == 2) then
          yhe4 = yhe4 + yk
       else
          ynuc = ynuc + yk
          a_num = a_num + net%a(k) * yk
          z_num = z_num + net%z(k) * yk
       endif
    enddo

    ! The NSE masses are bare-nuclear masses.  Add the rest mass of the
    ! balancing electrons so that the energy zero matches the convention
    ! used in statistic_compose() and in the DD2 <-> Helmholtz comparison.
    mexc = mexc + memev * ye_target

    if (ytot <= 0.d0) error stop "non-positive total ion abundance"
    abar = 1.d0 / ytot

    if (ynuc > 0.d0) then
       anuc = a_num / ynuc
       znuc = z_num / ynuc
    else
       anuc = 0.d0
       znuc = 0.d0
    endif

  end subroutine compose_moments


  subroutine fill_chemical_potentials(t, yq, q1, q6, mue, q3, q4, q5)
    real(8), intent(in) :: t(:), yq(:)
    real(8), intent(in) :: q1(:,:,:), q6(:,:,:), mue(:,:,:)
    real(8), intent(out) :: q3(:,:,:), q4(:,:,:), q5(:,:,:)

    integer :: inb, it, iyq
    integer :: nnb, nt, nyq
    real(8) :: dy
    real(8) :: mu_l, mu_q, mu_b, f_b
    real(8) :: mu_l_line(size(yq))

    nt  = size(t)
    nyq = size(yq)
    nnb = size(q1,3)

    dy = yq(2) - yq(1)
    if (dy <= 0.d0) error stop "yq grid must be strictly increasing"

    do iyq = 2, nyq
       if (abs((yq(iyq)-yq(iyq-1))/dy - 1.d0) > 1.d-10) &
            error stop "fill_chemical_potentials requires an equally spaced yq grid"
    enddo

    do inb = 1, nnb
       do it = 1, nt

          ! dF_b/dYq at fixed (nb,T), second order on the uniform Yq grid.
          mu_l_line(1) = mnmev * ( &
               -3.d0*(1.d0 + q6(it,1,inb)) &
               +4.d0*(1.d0 + q6(it,2,inb)) &
               -      (1.d0 + q6(it,3,inb)) ) / (2.d0*dy)

          do iyq = 2, nyq-1
             mu_l_line(iyq) = mnmev * &
                  (q6(it,iyq+1,inb) - q6(it,iyq-1,inb)) / (2.d0*dy)
          enddo

          mu_l_line(nyq) = mnmev * ( &
               +3.d0*(1.d0 + q6(it,nyq,inb)) &
               -4.d0*(1.d0 + q6(it,nyq-1,inb)) &
               +      (1.d0 + q6(it,nyq-2,inb)) ) / (2.d0*dy)

          do iyq = 1, nyq
             mu_l = mu_l_line(iyq)

             ! CompOSE: mu_e = mu_l - mu_q.
             mu_q = mu_l - mue(it,iyq,inb)

             ! F_b = -P/n_b + mu_b + Yq*mu_l.
             f_b = mnmev * (1.d0 + q6(it,iyq,inb))
             mu_b = f_b + q1(it,iyq,inb) - yq(iyq)*mu_l

             q3(it,iyq,inb) = mu_b / mnmev - 1.d0
             q4(it,iyq,inb) = mu_q / mnmev
             q5(it,iyq,inb) = mu_l / mnmev
          enddo

       enddo
    enddo

  end subroutine fill_chemical_potentials

end program make_compose_helmholtz
