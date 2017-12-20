! Basegrid.F90 : Register symmetries of the grid functions
!
!=============================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"


subroutine ScalarBase_symmetries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT ierr

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ScalarBase::phi1" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ScalarBase::phi2" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ScalarBase::Kphi1" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ScalarBase::Kphi2" )

end subroutine ScalarBase_symmetries
