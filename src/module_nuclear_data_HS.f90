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

    use iso_fortran_env, only: int32

    character(*), intent(in) :: fn

    integer :: iu, ios
    integer :: k, ia, iz
    integer(int32) :: record_size

    allocate(az_HS(nct_HS,2))
    allocate(mass_HS(nct_HS))
    allocate(bind_HS(nct_HS))

    open(newunit=iu, file=trim(fn), &
         status="old", action="read", &
         form="unformatted", access="stream", &
         convert="little_endian", &
         iostat=ios)

    if (ios /= 0) then
       write(*,*) "ERROR opening HS nuclear-data file: ", trim(fn)
       error stop
    endif

    ! Intel sequential-unformatted file:
    ! first 4 bytes = record length
    read(iu, pos=1, iostat=ios) record_size

    if (ios /= 0) then
       write(*,*) "ERROR reading HS record marker"
       error stop
    endif

    if (record_size /= 103284384_int32) then
       write(*,*) "ERROR: unexpected HS record size:", record_size
       error stop
    endif

    ! Payload begins immediately after the 4-byte record marker.
    ! Read only az, mass, bind; do not touch the huge comp array.
    read(iu, pos=5, iostat=ios) az_HS, mass_HS, bind_HS

    close(iu)

    if (ios /= 0) then
       write(*,*) "ERROR reading HS nuclear data: ", trim(fn)
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


  subroutine get_nuclear_data_HS(k, mass, bind, mexc)

    use const, only: mumev

    integer, intent(in) :: k
    real(8), intent(out) :: mass, bind, mexc

    if (k < 1 .or. k > nct_HS) then
       write(*,*) "ERROR: invalid HS nuclear-data index:", k
       error stop
    endif

    mass = mass_HS(k)
    bind = bind_HS(k)
    mexc = mass - dble(az_HS(k,1))*mumev

  end subroutine get_nuclear_data_HS

end module module_nuclear_data_HS
