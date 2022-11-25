
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

subroutine Scalar_zero_densities( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  integer                  istat
  logical                  use_jacobian

  call CCTK_IsFunctionAliased(istat, "MultiPatch_GetDomainSpecification")
  if (istat == 0) then
     use_jacobian = .false.
  else
     use_jacobian = .true.
  end if


  !--- Initialize fluxes grid functions (if statement needed to avoid nans in output) --------------
  if ( ( compute_fluxes == 1 ) .and. ( use_jacobian .eqv. .false. ) ) then
    rhoSF_gf      = 0.0d0
    jrSF_gf       = 0.0d0
    SrrSF_gf      = 0.0d0
  end if

end subroutine Scalar_zero_densities
