
#include "cctk.h"
#include "cctk_Arguments.h"

subroutine Scalar_zero_rhs( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS

  rhs_phi1    = 0
  rhs_phi2    = 0

  rhs_Kphi1   = 0
  rhs_Kphi2   = 0

end subroutine Scalar_zero_rhs
