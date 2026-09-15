module module_stat_weight_policy

  implicit none
  private

  integer, parameter, public :: STAT_WEIGHT_NONE     = 0
  integer, parameter, public :: STAT_WEIGHT_WINVNE   = 1
  integer, parameter, public :: STAT_WEIGHT_RAUSCHER = 2

  type, public :: stat_weight_policy_t
     integer :: primary  = STAT_WEIGHT_RAUSCHER
     integer :: fallback = STAT_WEIGHT_WINVNE
  end type stat_weight_policy_t

  type, public :: stat_weight_ref_t
     integer :: source = STAT_WEIGHT_NONE
     integer :: index  = 0
  end type stat_weight_ref_t
  
  public :: stat_weight_source_name
  public :: valid_stat_weight_policy
  public :: resolve_stat_weight_ref
  
contains

  function stat_weight_source_name(source) result(name)

    integer, intent(in) :: source
    character(len=16) :: name

    select case (source)

    case (STAT_WEIGHT_NONE)
       name = "none"

    case (STAT_WEIGHT_WINVNE)
       name = "winvne"

    case (STAT_WEIGHT_RAUSCHER)
       name = "rauscher"

    case default
       name = "unknown"

    end select

  end function stat_weight_source_name


  logical function valid_stat_weight_policy(policy)

    type(stat_weight_policy_t), intent(in) :: policy

    valid_stat_weight_policy = .false.

    select case (policy%primary)
    case (STAT_WEIGHT_WINVNE, STAT_WEIGHT_RAUSCHER)
       continue
    case default
       return
    end select

    select case (policy%fallback)
    case (STAT_WEIGHT_NONE, STAT_WEIGHT_WINVNE, STAT_WEIGHT_RAUSCHER)
       continue
    case default
       return
    end select

    if (policy%fallback == policy%primary) return

    valid_stat_weight_policy = .true.

  end function valid_stat_weight_policy

  subroutine resolve_stat_weight_ref(policy, iwinvne, irauscher, ref)

    type(stat_weight_policy_t), intent(in) :: policy
    integer, intent(in) :: iwinvne, irauscher
    type(stat_weight_ref_t), intent(out) :: ref

    ref%source = STAT_WEIGHT_NONE
    ref%index  = 0

    call try_stat_weight_source( &
         policy%primary, iwinvne, irauscher, ref)

    if (ref%source /= STAT_WEIGHT_NONE) return

    call try_stat_weight_source( &
         policy%fallback, iwinvne, irauscher, ref)

  end subroutine resolve_stat_weight_ref
  
  subroutine try_stat_weight_source(source, iwinvne, irauscher, ref)

    integer, intent(in) :: source
    integer, intent(in) :: iwinvne, irauscher
    type(stat_weight_ref_t), intent(inout) :: ref

    select case (source)

    case (STAT_WEIGHT_WINVNE)

       if (iwinvne > 0) then
          ref%source = STAT_WEIGHT_WINVNE
          ref%index  = iwinvne
       endif

    case (STAT_WEIGHT_RAUSCHER)

       if (irauscher > 0) then
          ref%source = STAT_WEIGHT_RAUSCHER
          ref%index  = irauscher
       endif

    case (STAT_WEIGHT_NONE)

       continue

    end select

  end subroutine try_stat_weight_source

end module module_stat_weight_policy
