/* ParamCheck.c : Check that the parameters provided make sense                  */
/* ============================================================================= */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void NPScalarsSF_ParamCheck(CCTK_ARGUMENTS){

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  for (int d = 0; d < cctk_dim; ++d)
      if (cctk_nghostzones[d] < (NP_order + 1) / 2)
          CCTK_ERROR ("This thorn requires at least (NP_order + 1) / 2 ghost zones");

}
