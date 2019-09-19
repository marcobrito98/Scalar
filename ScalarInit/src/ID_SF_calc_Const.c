/* intial data thorn: ID_SF_Const */
/*======================================================*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "ID_SF_utils.h"

/* -------------------------------------------------------------------*/
void ID_SF_Const(CCTK_ARGUMENTS);
void
ID_SF_Const (CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INFO("=== Begin ID_SF_Const initial data ===");

  /*=== define parameters ====================================*/
  /*==========================================================*/

  /*=== define grid length ===================================*/
  CCTK_INT imin[3], imax[3];
  for (int d = 0; d < 3; ++ d)
  {
    imin[d] = 0;
    imax[d] = cctk_lsh[d];
  }

  /*==========================================================*/

/*=== loops over full grid ===================================*/
/*------------------------------------------------------------*/
//#pragma omp parallel for
  for (int k = imin[2]; k < imax[2]; ++k)
  {
   for (int j = imin[1]; j < imax[1]; ++j)
   {
    for (int i = imin[0]; i < imax[0]; ++i)
    {

     const int ind = CCTK_GFINDEX3D (cctkGH, i, j, k);

     /*=== initialize grid functions as zero ==================*/
     phi1[ind]  = 0.0;
     phi2[ind]  = 0.0;
     Kphi1[ind] = 0.0;
     Kphi2[ind] = 0.0;
     /*========================================================*/

     /*=== initialize local functions as zero =================*/
     // scalar field momentum
     CCTK_REAL psit_re, psit_im;
     psit_re = 0.0;
     psit_im = 0.0;
     /*========================================================*/

     /*=== define radius and angles ===========================*/
     // positions OLD with KerrQuasiIsotropic
     /*CCTK_REAL xp[3];
     xp[0] = x[ind] - pos_plus[0];
     xp[1] = y[ind] - pos_plus[1];
     xp[2] = z[ind] - pos_plus[2];*/

     // positions NEW
     CCTK_REAL xp[3];
     xp[0] = x[ind];
     xp[1] = y[ind];
     xp[2] = z[ind];

     // coordinate radius and polar radial coordinate
     CCTK_REAL rr, rr2;
     rr2 = xp[0] * xp[0] + xp[1] * xp[1] + xp[2] * xp[2];
     if( rr2 < eps_r2 ) rr2 = eps_r2;
     rr  = sqrt( rr2 );

     CCTK_REAL rho, rho2;
     rho2 = xp[0]*xp[0] + xp[1]*xp[1];
     if( rho2 < eps_r2 ) rho2 = eps_r2;
     rho  = sqrt( rho2 );


     // angles
     CCTK_REAL ctheta, ctheta2, cphi, cphi2;
     CCTK_REAL stheta, stheta2, sphi, sphi2;

     ctheta  = xp[2] / rr;
     ctheta2 = ctheta * ctheta;
     cphi    = xp[0] / rho;
     cphi2   = cphi * cphi;

     stheta  = rho / rr;
     stheta2 = stheta * stheta;
     sphi    = xp[1] / rho;
     sphi2   = sphi * sphi;

     CCTK_REAL c2phi, s2phi;
     c2phi   = cphi2 - sphi2;
     s2phi   = 2.0 * cphi * sphi;
     /*========================================================*/

     /*=== initialize spherical harmonics ( in cartesian coordinates ) ===*/

     /*--------------------------------------------------------*/
	if( l0SF == 0 )
	// l0SF = m0SF = 0
	{
	  psit_re = 1.0 / sqrt( 4.0*Pi );
	  psit_im = 0.0;
	}
	/*-----------------------------------------------------*/
	else if( l0SF == 1 && m0SF == -1 )
	// l0SF = 1, m0SF = -1
	{
	  psit_re =   sqrt( 3.0 / ( 8.0*Pi ) ) * cphi * stheta;
	  psit_im = - sqrt( 3.0 / ( 8.0*Pi ) ) * sphi * stheta;
	}
	else if( l0SF == 1 && m0SF == 0)
	// l0SF = 1, m0SF = 0
	{
	  psit_re = sqrt( 3.0 / ( 4.0*Pi ) ) * ctheta;
	  psit_im = 0.0;
	}
	else if( l0SF == 1 && m0SF == 1)
	// l0SF = 1, m0SF = 1
	{
	  psit_re = - sqrt( 3.0 / ( 8.0*Pi ) ) * cphi * stheta;
	  psit_im = - sqrt( 3.0 / ( 8.0*Pi ) ) * sphi * stheta;
	}
	/*-----------------------------------------------------*/
	else if( l0SF == 2 && m0SF == -2 )
	// l0SF = 2, m0SF = -2
	{ 
	  psit_re =   sqrt( 15.0 / ( 32.0*Pi ) ) * c2phi * stheta2;
	  psit_im = - sqrt( 15.0 / ( 32.0*Pi ) ) * s2phi * stheta2;
	}
	else if( l0SF == 2 && m0SF == -1 )
	// l0SF = 2, m0SF = -1
	{ 
	  psit_re =   sqrt( 15.0 / (  8.0*Pi ) ) * cphi * ctheta * stheta;
	  psit_im = - sqrt( 15.0 / (  8.0*Pi ) ) * sphi * ctheta * stheta;
	}
	else if( l0SF == 2 && m0SF == 0 )
	// l0SF = 2, m0SF = 0
	{ 
	  psit_re =   sqrt(  5.0 / ( 16.0*Pi ) ) * ( 3.0 * ctheta2 - 1.0 );
	  psit_im = 0.0;
	}
	else if( l0SF == 2 && m0SF == 1 )
	// l0SF = 2, m0SF = 1
	{ 
	  psit_re = - sqrt( 15.0 / (  8.0*Pi ) ) * cphi * ctheta * stheta;
	  psit_im = - sqrt( 15.0 / (  8.0*Pi ) ) * sphi * ctheta * stheta;
	}
	else if( l0SF == 2 && m0SF == 2 )
	// l0SF = 2, m0SF = 2
	{
	  psit_re =   sqrt( 15.0 / ( 32.0*Pi ) ) * c2phi * stheta2;
	  psit_im =   sqrt( 15.0 / ( 32.0*Pi ) ) * s2phi * stheta2;
	}
	/*-----------------------------------------------------*/
     	else
     	CCTK_WARN (0, "invalid multipole for scalar field initial data");
     /*========================================================*/

     /*=== calc momentum Kphi =================================*/

     psit_re = ampSF * psit_re; 
     psit_im = ampSF * psit_im;

     /*========================================================*/

     /*=== write to grid functions ============================*/

     phi1[ind]  = psit_re;
     phi2[ind]  = psit_im;

     Kphi1[ind] = 0;
     Kphi2[ind] = 0;

     /*========================================================*/

    } /* for i */
   }  /* for j */
  }   /* for k */
/*=== end of loops over grid =================================*/
/*============================================================*/

  CCTK_INFO("=== End ID_SF_Const initial data ===");

}
/* -------------------------------------------------------------------*/
