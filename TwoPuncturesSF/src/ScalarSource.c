/* TwoPuncturesSF:  File  "ScalarSource.c"*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk_Parameters.h"
#include "TP_utilities.h"
#include "TwoPuncturesSF.h"


void SF_Gaussian(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL *phi2, CCTK_REAL *dphi2)
{
  DECLARE_CCTK_PARAMETERS

  /*=== initialize local functions as zero =================*/
  // scalar field momentum
  CCTK_REAL psit_re, psit_im;
  psit_re = 0.0;
  psit_im = 0.0;
  /*========================================================*/

  /*=== define radius and angles ===========================*/
  // positions
  CCTK_REAL xp[3];
  xp[0] = x - par_b; 
  xp[1] = y;
  xp[2] = z;

  // coordinate radius and polar radial coordinate
  CCTK_REAL rr, rr2;
  rr2 = xp[0] * xp[0] + xp[1] * xp[1] + xp[2] * xp[2];
  if( rr2 < TP_Tiny ) rr2 = TP_Tiny;
  rr  = sqrt( rr2 );

  /*========================================================*/


  if( CCTK_Equals(scalar_GaussProfile, "single_mode"))
  // single mode data
  {
    if( l0SF == 0 )
    // l0SF = m0SF = 0
    {
    psit_re = 1.0 / sqrt( 4.0*Pi );
    psit_im = 0.0;
    }
    /*-----------------------------------------------------*/
    else
    CCTK_WARN (0, "invalid multipole for scalar field initial data");
  }
  /*--------------------------------------------------------*/
  else
  CCTK_WARN (0, "invalid scalar field initial data");
  /*========================================================*/

  /*=== calc scalar field phi =================================*/

  *phi2 = pow( ampSF * psit_re, 2) * exp( - 2.0 * ( rr - r0SF )*( rr - r0SF ) / ( widthSF*widthSF ) ); 

  *dphi2 = pow( ampSF * psit_re * 2.0 * ( rr - r0SF ) / ( widthSF*widthSF ), 2)
      * exp( - 2.0 * ( rr - r0SF )*( rr - r0SF ) / ( widthSF*widthSF ) );

}

CCTK_REAL NonLinSrcSF (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL psi)
{
  DECLARE_CCTK_PARAMETERS;
  CCTK_REAL psi2, psi4, psi5;
  CCTK_REAL phi2, dphi2, term1, term2;

  psi2 = psi * psi;
  psi4 = psi2 * psi2;
  psi5 = psi * psi4;

  SF_Gaussian(x, y, z, &phi2, &dphi2);
  term1 = mu * mu * phi2 * psi5;
  term2 = dphi2 * psi;

  return Pi * (term1 + term2);
}

CCTK_REAL LinSrcSF (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL psi)
{
  DECLARE_CCTK_PARAMETERS;
  CCTK_REAL psi2, psi4;
  CCTK_REAL phi2, dphi2, term1, term2;

  psi2 = psi * psi;
  psi4 = psi2 * psi2;

  SF_Gaussian(x, y, z, &phi2, &dphi2);
  term1 = 5 * mu * mu * phi2 * psi4;
  term2 = dphi2;

  return Pi * (term1 + term2);
}

/*-----------------------------------------------------------*/
