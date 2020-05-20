#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

void Scalar_Init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // running with multipatch?
  if (!CCTK_IsFunctionAliased("MultiPatch_GetDomainSpecification") && z_is_radial) {
    CCTK_ERROR("Without multipatch, please set z_is_radial to no.");
    return;
  }

  if (CCTK_Equals(outer_bound, "reflecting") && !z_is_radial)
    CCTK_ERROR("this outer_bound only runs with Llama.");

  if (z_is_radial && !CCTK_IsThornActive("Interpolate2"))
    CCTK_ERROR("Please activate thorn Interpolate2 to use this boundary condition.");

    return;
}
