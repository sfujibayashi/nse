module module_nse_species_policy

  implicit none
  private

  integer, parameter, public :: NSE_SPECIES_LEGACY     = 1
  integer, parameter, public :: NSE_SPECIES_ALL_WINVNE = 2

  type, public :: nse_species_policy_t
     integer :: mode = NSE_SPECIES_LEGACY
  end type nse_species_policy_t

  public :: keep_nse_species
  public :: valid_nse_species_policy

contains

  logical function keep_nse_species(policy, iz, has_rauscher)

    type(nse_species_policy_t), intent(in) :: policy
    integer, intent(in) :: iz
    logical, intent(in) :: has_rauscher

    select case (policy%mode)

    case (NSE_SPECIES_LEGACY)

       keep_nse_species = has_rauscher .or. (iz < 87)

    case (NSE_SPECIES_ALL_WINVNE)

       keep_nse_species = .true.

    case default

       keep_nse_species = .false.

    end select

  end function keep_nse_species


  logical function valid_nse_species_policy(policy)

    type(nse_species_policy_t), intent(in) :: policy

    select case (policy%mode)

    case (NSE_SPECIES_LEGACY, NSE_SPECIES_ALL_WINVNE)
       valid_nse_species_policy = .true.

    case default
       valid_nse_species_policy = .false.

    end select

  end function valid_nse_species_policy

end module module_nse_species_policy
