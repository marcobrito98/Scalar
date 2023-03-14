! NPScalars_SF
! NPboundaries.F90: define symmetry boundaries for NP scalars
!
!=============================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine NPSFboundaries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT ierr

  CCTK_INT, parameter :: one = 1
  CCTK_INT width

  if (calculate_NP_every .le. 0) then
     return
  end if

  if (MOD(cctk_iteration, calculate_NP_every) .ne. 0 ) then
     return
  endif

  ! let's just for simplicity say that
  width = 1

  ! since these are just used for analysis, register all BCs as "flat"

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, width, -one,      &
       "NPScalars_SF::NPPsi4R_group", "flat")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for NPScalars::NPPsi4R_group!")

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, width, -one,      &
       "NPScalars_SF::NPPsi4I_group", "flat")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for NPScalars_SF::NPPsi4I_group!")

!  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, width, -one,      &
!       "NPScalars_SF::ji", "flat")
!  if (ierr < 0)                                                            &
!       call CCTK_WARN(0, "Failed to register BC for NPScalars_SF::ji!")

!  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, width, -one,      &
!       "NPScalars_SF::Sij", "flat")
!  if (ierr < 0)                                                            &
!       call CCTK_WARN(0, "Failed to register BC for NPScalars_SF::Sij!")


end subroutine NPSFboundaries
