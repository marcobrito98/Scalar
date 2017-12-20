#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

subroutine Scalar_InitSymBound( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT ierr

  ! there is no need to register the variables phi and Kphi themselves since
  ! this is done in the 'Base' thorn that defines them

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ScalarEvolve::rhs_phi1" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ScalarEvolve::rhs_phi2" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ScalarEvolve::rhs_Kphi1" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ScalarEvolve::rhs_Kphi2" )

end subroutine Scalar_InitSymBound
