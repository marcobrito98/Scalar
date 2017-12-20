
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void Scalar_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0, group, rhs;

  // register evolution and rhs gridfunction groups with MoL

  /* ADM metric and extrinsic curvature */
  group = CCTK_GroupIndex("ADMBase::lapse");
  ierr += MoLRegisterSaveAndRestoreGroup(group);
  group = CCTK_GroupIndex("ADMBase::shift");
  ierr += MoLRegisterSaveAndRestoreGroup(group);
  group = CCTK_GroupIndex("ADMBase::metric");
  ierr += MoLRegisterSaveAndRestoreGroup(group);
  group = CCTK_GroupIndex("ADMBase::curv");
  ierr += MoLRegisterSaveAndRestoreGroup(group);

  /* phi and rhs_phi */
  group = CCTK_GroupIndex("ScalarBase::phi");
  rhs   = CCTK_GroupIndex("ScalarEvolve::rhs_phi");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* Kphi and rhs_Kphi */
  group = CCTK_GroupIndex("ScalarBase::Kphi");
  rhs   = CCTK_GroupIndex("ScalarEvolve::rhs_Kphi");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  if (ierr) CCTK_ERROR("Problems registering with MoL");

}
