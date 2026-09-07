module module_ptf_rauscher

  implicit none

  integer,parameter,public :: nt_rauscher = 72

  real(8),parameter :: t9_rauscher(nt_rauscher) = (/ &
       0.1d0,  0.15d0, 0.2d0,  0.3d0,  0.4d0,  0.5d0,  &
       0.6d0,  0.7d0,  0.8d0,  0.9d0,  1.0d0,  1.5d0,  &
       2.0d0,  2.5d0,  3.0d0,  3.5d0,  4.0d0,  4.5d0,  &
       5.0d0,  6.0d0,  7.0d0,  8.0d0,  9.0d0, 10.0d0,  &
      12.0d0, 14.0d0, 16.0d0, 18.0d0, 20.0d0, 22.0d0,  &
      24.0d0, 26.0d0, 28.0d0, 30.0d0, 35.0d0, 40.0d0,  &
      45.0d0, 50.0d0, 55.0d0, 60.0d0, 65.0d0, 70.0d0,  &
      75.0d0, 80.0d0, 85.0d0, 90.0d0, 95.0d0,100.0d0,  &
     105.0d0,110.0d0,115.0d0,120.0d0,125.0d0,130.0d0,  &
     135.0d0,140.0d0,145.0d0,150.0d0,155.0d0,160.0d0,  &
     165.0d0,170.0d0,175.0d0,180.0d0,190.0d0,200.0d0,  &
     210.0d0,220.0d0,230.0d0,240.0d0,250.0d0,275.0d0 /)

  integer :: nct_rauscher = 0

  integer,allocatable :: z_rauscher(:)
  integer,allocatable :: a_rauscher(:)

  real(8),allocatable :: spin_rauscher(:)
  real(8),allocatable :: ptf_rauscher(:,:)

  public :: init_ptf_rauscher

contains
  
  subroutine init_ptf_rauscher(fn)

    character(*),intent(in) :: fn

    integer :: iu
    integer :: ios
    integer :: i, k
    integer :: iz, ia
    real(8) :: j0
    real(8) :: pf(nt_rauscher)

    character(1024) :: line

    ! ---------------------------------------------------------
    ! First pass: count nuclei
    ! ---------------------------------------------------------

    open(newunit=iu,file=fn,status="old",action="read")

    ! skip 23-line header
    do i=1,23
       read(iu,'(A)',iostat=ios) line
       if (ios /= 0) then
          write(*,*) "ERROR: failed while reading Rauscher header"
          stop
       endif
    enddo

    nct_rauscher = 0

    do
       read(iu,'(A)',iostat=ios) line
       if (ios /= 0) exit

       if (len_trim(line) == 0) cycle

       nct_rauscher = nct_rauscher + 1
    enddo

    close(iu)

    write(*,*) "# of Rauscher species:",nct_rauscher

    ! ---------------------------------------------------------
    ! Allocate
    ! ---------------------------------------------------------

    allocate(z_rauscher(nct_rauscher))
    allocate(a_rauscher(nct_rauscher))
    allocate(spin_rauscher(nct_rauscher))
    allocate(ptf_rauscher(nt_rauscher, nct_rauscher))

    ! ---------------------------------------------------------
    ! Second pass: read data
    ! ---------------------------------------------------------

    open(newunit=iu,file=fn,status="old",action="read")

    do i=1,23
       read(iu,'(A)') line
    enddo

    k = 0

    do
       read(iu,'(A)',iostat=ios) line
       if (ios /= 0) exit

       if (len_trim(line) == 0) cycle

       read(line,*,iostat=ios) iz, ia, j0, pf(:)

       if (ios /= 0) then
          write(*,*) "ERROR reading Rauscher record:"
          write(*,'(A)') trim(line)
          stop
       endif

       k = k + 1

       z_rauscher(k) = iz
       a_rauscher(k) = ia
       spin_rauscher(k) = j0
       ptf_rauscher(:,k) = pf(:)

    enddo

    close(iu)

    if (k /= nct_rauscher) then
       write(*,*) "ERROR: inconsistent species count:", &
            k,nct_rauscher
       stop
    endif

  end subroutine init_ptf_rauscher

end module module_ptf_rauscher
