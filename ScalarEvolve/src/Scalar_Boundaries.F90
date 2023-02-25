#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine Scalar_Boundaries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT ierr

  CCTK_INT, parameter :: one = 1

  ! The outgoing (radiative) boundary conditions are being handled from calls to
  ! the NewRad infraestructure. Here we register all BCs as 'none', which
  ! enforces all the symmetry BCs.

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,      &
       "ScalarBase::phi", "none")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for ScalarBase::phi!")

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,      &
       "ScalarBase::Kphi", "none")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for ScalarBase::Kphi!")

end subroutine Scalar_Boundaries

!=============================================================================
!
subroutine Scalar_calc_Tmunu_bdry( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT ierr
  CCTK_INT, parameter :: one = 1
  CCTK_INT bndsize

  bndsize = 3

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, bndsize, -one, &
       "TmunuBase::stress_energy_scalar", "flat")
  if (ierr < 0)                                                           &
       call CCTK_ERROR("Failed to register BC for TmunuBase::stress_energy_scalar!")

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, bndsize, -one, &
       "TmunuBase::stress_energy_vector", "flat")
  if (ierr < 0)                                                           &
       call CCTK_ERROR("Failed to register BC for TmunuBase::stress_energy_vector!")

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, bndsize, -one, &
       "TmunuBase::stress_energy_tensor", "flat")
  if (ierr < 0)                                                           &
       call CCTK_ERROR("Failed to register BC for TmunuBase::stress_energy_tensor!")

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, bndsize, -one, &
       "ScalarEvolve::densities_SF", "flat")
  if (ierr < 0)                                                           &
       call CCTK_ERROR("Failed to register BC for ScalarEvolve::densities_SF!")

end subroutine Scalar_calc_Tmunu_bdry
