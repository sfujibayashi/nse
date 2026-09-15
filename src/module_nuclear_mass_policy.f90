module module_nuclear_mass_policy

  implicit none
  private

  integer, parameter, public :: NUCLEAR_MASS_NONE   = 0
  integer, parameter, public :: NUCLEAR_MASS_WINVNE = 1
  integer, parameter, public :: NUCLEAR_MASS_HS     = 2

  type, public :: nuclear_mass_policy_t
     integer :: primary  = NUCLEAR_MASS_WINVNE
     integer :: fallback = NUCLEAR_MASS_NONE
  end type nuclear_mass_policy_t

  type, public :: nuclear_mass_ref_t
     integer :: source = NUCLEAR_MASS_NONE
     integer :: index  = 0
  end type nuclear_mass_ref_t

  public :: nuclear_mass_source_name
  public :: valid_nuclear_mass_policy
  public :: resolve_nuclear_mass_ref

contains

  function nuclear_mass_source_name(source) result(name)

    integer, intent(in) :: source
    character(len=16) :: name

    select case (source)

    case (NUCLEAR_MASS_NONE)
       name = "none"

    case (NUCLEAR_MASS_WINVNE)
       name = "winvne"

    case (NUCLEAR_MASS_HS)
       name = "HS"

    case default
       name = "unknown"

    end select

  end function nuclear_mass_source_name


  logical function valid_nuclear_mass_policy(policy)

    type(nuclear_mass_policy_t), intent(in) :: policy

    valid_nuclear_mass_policy = .false.

    select case (policy%primary)

    case (NUCLEAR_MASS_WINVNE, NUCLEAR_MASS_HS)
       continue

    case default
       return

    end select

    select case (policy%fallback)

    case (NUCLEAR_MASS_NONE, NUCLEAR_MASS_WINVNE, NUCLEAR_MASS_HS)
       continue

    case default
       return

    end select

    if (policy%fallback == policy%primary) return

    valid_nuclear_mass_policy = .true.

  end function valid_nuclear_mass_policy


  subroutine resolve_nuclear_mass_ref( &
       policy, iwinvne, ihs, ref)

    type(nuclear_mass_policy_t), intent(in) :: policy
    integer, intent(in) :: iwinvne, ihs
    type(nuclear_mass_ref_t), intent(out) :: ref

    ref%source = NUCLEAR_MASS_NONE
    ref%index  = 0

    call try_nuclear_mass_source( &
         policy%primary, iwinvne, ihs, ref)

    if (ref%source /= NUCLEAR_MASS_NONE) return

    call try_nuclear_mass_source( &
         policy%fallback, iwinvne, ihs, ref)

  end subroutine resolve_nuclear_mass_ref


  subroutine try_nuclear_mass_source( &
       source, iwinvne, ihs, ref)

    integer, intent(in) :: source
    integer, intent(in) :: iwinvne, ihs
    type(nuclear_mass_ref_t), intent(inout) :: ref

    select case (source)

    case (NUCLEAR_MASS_WINVNE)

       if (iwinvne > 0) then
          ref%source = NUCLEAR_MASS_WINVNE
          ref%index  = iwinvne
       endif

    case (NUCLEAR_MASS_HS)

       if (ihs > 0) then
          ref%source = NUCLEAR_MASS_HS
          ref%index  = ihs
       endif

    case (NUCLEAR_MASS_NONE)

       continue

    end select

  end subroutine try_nuclear_mass_source

end module module_nuclear_mass_policy
