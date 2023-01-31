
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

/* -------------------------------------------------------------------*/
void
ID_SF_BS_ParamCheck(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if( m_plus != 1.0 )
    {
        CCTK_WARN( 0, "BH mass value not valid. Bound state initial data works only with m_plus = 1 ");
    }
}
