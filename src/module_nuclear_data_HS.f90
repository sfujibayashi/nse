module module_nuclear_data_HS

  implicit none
  private

  integer, parameter, public :: nct_HS = 8140

  integer, allocatable :: az_HS(:,:)
  real(8), allocatable :: mass_HS(:)
  real(8), allocatable :: bind_HS(:)

  integer, allocatable :: jnuc_HS(:,:)

  integer :: na_HS = 0
  integer :: nz_HS = 0

  public :: init_nuclear_data_HS
  public :: find_HS_index
  public :: get_nuclear_data_HS

contains

  subroutine init_nuclear_data_HS(fn)

    character(*), intent(in) :: fn

    integer :: iu, ios
    integer :: k, ia, iz

    allocate(az_HS(nct_HS,2))
    allocate(mass_HS(nct_HS))
    allocate(bind_HS(nct_HS))

    open(newunit=iu, file=trim(fn), &
         status="old", action="read", &
         form="unformatted", access="sequential", &
         iostat=ios)

    if (ios /= 0) then
       write(*,*) "ERROR opening HS nuclear-data file: ", trim(fn)
       error stop
    endif

    ! The original record contains
    !
    !   az, mass, bind, comp
    !
    ! Only the first three objects are needed here.
    read(iu, iostat=ios) az_HS, mass_HS, bind_HS

    close(iu)

    if (ios /= 0) then
       write(*,*) "ERROR reading HS nuclear-data file: ", trim(fn)
       error stop
    endif

    na_HS = maxval(az_HS(:,1))
    nz_HS = maxval(az_HS(:,2))

    allocate(jnuc_HS(1:na_HS,0:nz_HS))
    jnuc_HS(:,:) = 0

    do k = 1, nct_HS

       ia = az_HS(k,1)
       iz = az_HS(k,2)

       if (ia < 1 .or. iz < 0) cycle

       if (jnuc_HS(ia,iz) /= 0) then
          write(*,*) "ERROR: duplicate HS nucleus:", ia, iz
          error stop
       endif

       jnuc_HS(ia,iz) = k

    enddo

    write(*,*) "# of HS species:", nct_HS
    write(*,*) "HS Z, A max:", nz_HS, na_HS

  end subroutine init_nuclear_data_HS


  integer function find_HS_index(ia, iz) result(idx)

    integer, intent(in) :: ia, iz

    idx = 0

    if (.not. allocated(jnuc_HS)) return

    if (ia < lbound(jnuc_HS,1) .or. &
         ia > ubound(jnuc_HS,1)) return

    if (iz < lbound(jnuc_HS,2) .or. &
         iz > ubound(jnuc_HS,2)) return

    idx = jnuc_HS(ia,iz)

  end function find_HS_index


  subroutine get_nuclear_data_HS(k, mass, bind)

    integer, intent(in) :: k
    real(8), intent(out) :: mass, bind

    if (k < 1 .or. k > nct_HS) then
       write(*,*) "ERROR: invalid HS nuclear-data index:", k
       error stop
    endif

    mass = mass_HS(k)
    bind = bind_HS(k)

  end subroutine get_nuclear_data_HS

end module module_nuclear_data_HS
