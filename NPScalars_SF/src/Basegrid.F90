! NPScalars_SF
! Basegrid.F90 : Register symmetries
!
!=============================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

subroutine NPSF_symmetries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT  ierr

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "NPScalars_SF::psi4re" )
  call SetCartSymVN( ierr, cctkGH, (/-1,-1,-1/), "NPScalars_SF::psi4im" )

end subroutine NPSF_symmetries
