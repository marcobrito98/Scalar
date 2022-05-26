#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "TP_utilities.h"
#include "TwoPunctures_BBHSF.h"
#include <assert.h>

/* -------------------------------------------------------------------*/
void
SF_Initialize(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL *phi_re, CCTK_REAL *phi_im)
{

  DECLARE_CCTK_PARAMETERS;

     /*=== define radius and angles ===========================*/
     // positions
     CCTK_REAL xp[3];
     xp[0] = x;
     //xp[0] = x - par_b;
     xp[1] = y;
     xp[2] = z;

#define EXTEND(M,r) \
          ( M * (3./8 * pow(r, 4) / pow(TP_Extend_Radius, 5) - \
                 5./4 * pow(r, 2) / pow(TP_Extend_Radius, 3) + \
                 15./8 / TP_Extend_Radius))

     // coordinate radius rr and polar radial coordinate rho
     CCTK_REAL rr, rr2, irr, irr2;
     rr = sqrt(xp[0] * xp[0] + xp[1] * xp[1] + xp[2] * xp[2]);
     rr = pow (pow (rr, 4) + pow (TP_epsilon, 4), 0.25);
     if( rr < TP_Tiny ) rr = TP_Tiny;
     irr = 1.0/rr;
     if (abs(rr) < TP_Extend_Radius) {
           irr = EXTEND(1., rr);
     }
     rr2 = rr * rr;
     irr2  = irr  * irr;
 
     CCTK_REAL rho, rho2, irho, irho2;
     rho = sqrt( xp[0]*xp[0] + xp[1]*xp[1] );
     rho = pow (pow (rho, 4) + pow (TP_epsilon, 4), 0.25);
     if( rho < TP_Tiny ) rho = TP_Tiny;
     irho = 1.0/rho;
     if (abs(rho) < TP_Extend_Radius) {
           irho = EXTEND(1., rho);
     }
     rho2   = rho * rho;
     irho2  = irho * irho;

     CCTK_REAL deltR =  rr - r0SF;

     // angles
     CCTK_REAL ctheta, ctheta2, cphi, cphi2;
     CCTK_REAL stheta, stheta2, sphi, sphi2;

     ctheta  = xp[2] * irr;
     ctheta2 = ctheta * ctheta;
     cphi    = xp[0] * irho;
     cphi2   = cphi * cphi;

     stheta  = rho * irr;
     stheta2 = stheta * stheta;
     sphi    = xp[1] * irho;
     sphi2   = sphi * sphi;

     CCTK_REAL c2phi, s2phi;
     c2phi   = cphi2 - sphi2;
     s2phi   = 2.0 * cphi * sphi;

     /*========================================================*/
        CCTK_INT monopole_statement, dipole_statement;

        /* gaussian l=m=0 scalar field  */
        monopole_statement = (CCTK_EQUALS(initial_scalar, "ID_SF_Gaussian") && CCTK_EQUALS(scalar_GaussProfile, "single_mode")) && ( l0SF == 0 ) && ( m0SF == 0 );
        /* gaussian l=m=1 scalar field  */
        dipole_statement = (CCTK_EQUALS(initial_scalar, "ID_SF_Gaussian") && CCTK_EQUALS(scalar_GaussProfile, "single_mode")) && ( l0SF == 1 ) && ( m0SF == 1 );

        if( monopole_statement )
        {
                *phi_re = (1. / sqrt(4*Pi)) * ampSF * exp( - deltR * deltR / ( widthSF * widthSF ) );
                *phi_im = 0.;
        }

	else if( dipole_statement )
        {
                *phi_re = - ampSF / 2. * sqrt(1.5 / Pi) * stheta * cphi * exp( - deltR*deltR / ( widthSF*widthSF ) );
                *phi_im = - ampSF / 2. * sqrt(1.5 / Pi) * stheta * sphi * exp( - deltR*deltR / ( widthSF*widthSF ) );
        }

	else
	{
		CCTK_WARN (0, "TwoPunctures_BBHSF: invalid choice for scalar field initial profile");	
	}

}


void conf_flat_analytic_SF_source_term (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL phi_re, CCTK_REAL phi_im, CCTK_REAL *sourceSF1, CCTK_REAL *sourceSF2, CCTK_REAL Dphi_re_dx[3], CCTK_REAL Dphi_im_dx[3])
{
  DECLARE_CCTK_PARAMETERS;

  /*=== define radius and angles ===========================*/
  // positions
  CCTK_REAL xp[3];
  xp[0] = x;
  //xp[0] = x - par_b;
  xp[1] = y;
  xp[2] = z;

#define EXTEND(M,r) \
          ( M * (3./8 * pow(r, 4) / pow(TP_Extend_Radius, 5) - \
                 5./4 * pow(r, 2) / pow(TP_Extend_Radius, 3) + \
                 15./8 / TP_Extend_Radius))
 
  // coordinate radius rr and polar radial coordinate rho
  CCTK_REAL rr, rr2, irr, irr2;
  rr = sqrt(xp[0] * xp[0] + xp[1] * xp[1] + xp[2] * xp[2]);
  rr = pow (pow (rr, 4) + pow (TP_epsilon, 4), 0.25);
  if( rr < TP_Tiny ) rr = TP_Tiny;
  irr = 1.0/rr;
  if (abs(rr) < TP_Extend_Radius) {
     irr = EXTEND(1., rr);
  }
  rr2 = rr * rr;
  irr2  = irr  * irr;

  CCTK_REAL rho, rho2, irho, irho2;
  rho = sqrt(xp[0]*xp[0] + xp[1]*xp[1]);
  rho = pow (pow (rho, 4) + pow (TP_epsilon, 4), 0.25);
  if( rho < TP_Tiny ) rho = TP_Tiny;
  irho = 1.0/rho;
  if (abs(rho) < TP_Extend_Radius) {
     irho = EXTEND(1., rho);
  }
  rho2   = rho * rho;
  irho2  = irho * irho;

  CCTK_REAL deltR =  rr - r0SF;

  // angles
  CCTK_REAL ctheta, ctheta2, cphi, cphi2;
  CCTK_REAL stheta, stheta2, sphi, sphi2;

  ctheta  = xp[2] * irr;
  ctheta2 = ctheta * ctheta;
  cphi    = xp[0] * irho;
  cphi2   = cphi * cphi;

  stheta  = rho * irr;
  stheta2 = stheta * stheta;
  sphi    = xp[1] * irho;
  sphi2   = sphi * sphi;

  CCTK_REAL c2phi, s2phi;
  c2phi   = cphi2 - sphi2;
  s2phi   = 2.0 * cphi * sphi;

  /*--- spherical coordinates to cartesian coordinates derivatives ---*/
  CCTK_REAL dRdx, dRdy, dRdz, dTdx, dTdy, dTdz, dPdx, dPdy;
  dRdx =   xp[0] * irr;
  dRdy =   xp[1] * irr;
  dRdz =   xp[2] * irr;

  dTdx =   xp[0] * xp[2] * irho * irr2;
  dTdy =   xp[1] * xp[2] * irho * irr2;
  dTdz = - rho * irr2;

  dPdx = - xp[1] * irho2;
  dPdy =   xp[0] * irho2;

  /*========================================================*/

  // SF source term
  CCTK_INT monopole_statement, dipole_statement;
        
  /* gaussian l=m=0 scalar field choice*/
  monopole_statement = (CCTK_EQUALS(initial_scalar, "ID_SF_Gaussian") && CCTK_EQUALS(scalar_GaussProfile, "single_mode")) && ( l0SF == 0 ) && ( m0SF == 0 );
  /* gaussian l=m=1 scalar field choice*/
  dipole_statement = (CCTK_EQUALS(initial_scalar, "ID_SF_Gaussian") && CCTK_EQUALS(scalar_GaussProfile, "single_mode")) && ( l0SF == 1 ) && ( m0SF == 1 );

  /*--- scalar field ---*/
  CCTK_REAL psi_re, psi_im;
  CCTK_INT  i;

  psi_re = phi_re;
  psi_im = phi_im;

       if ( monopole_statement ) {

        	CCTK_REAL DpsiRdR;
        	CCTK_REAL DpsiIdR;

        	DpsiRdR = - 2. * deltR * psi_re / ( widthSF * widthSF );
        	DpsiIdR = 0.;

		Dphi_re_dx[0] = DpsiRdR * dRdx;
		Dphi_re_dx[1] = DpsiRdR * dRdy;
		Dphi_re_dx[2] = DpsiRdR * dRdz;

		Dphi_im_dx[0] = DpsiIdR * dRdx;
		Dphi_im_dx[1] = DpsiIdR * dRdy;
		Dphi_im_dx[2] = DpsiIdR * dRdz;

	        for( i=0; i<3; i++){	
                *sourceSF1 += Dphi_re_dx[i] * Dphi_re_dx[i] + Dphi_im_dx[i] * Dphi_im_dx[i];
		}

                *sourceSF2 = mu * mu * ( psi_re * psi_re + psi_im * psi_im );
	}
	
	else if ( dipole_statement ) {

                CCTK_REAL DpsiRdR, DpsiRdT, DpsiRdP;
                CCTK_REAL DpsiIdR, DpsiIdT, DpsiIdP;

                DpsiRdR = - 2. * deltR * psi_re / ( widthSF * widthSF );
                DpsiRdT = psi_re * ctheta / stheta;
                DpsiRdP = - psi_re * sphi / cphi;
                
		DpsiIdR = - 2. * deltR * psi_im / ( widthSF * widthSF );
		DpsiIdT = psi_im * ctheta / stheta;
		DpsiIdP = psi_im * cphi / sphi;

		Dphi_re_dx[0] = DpsiRdR * dRdx + DpsiRdT * dTdx + DpsiRdP * dPdx;
        	Dphi_re_dx[1] = DpsiRdR * dRdy + DpsiRdT * dTdy + DpsiRdP * dPdy;
        	Dphi_re_dx[2] = DpsiRdR * dRdz + DpsiRdT * dTdz ;

        	Dphi_im_dx[0] = DpsiIdR * dRdx + DpsiIdT * dTdx + DpsiIdP * dPdx;
        	Dphi_im_dx[1] = DpsiIdR * dRdy + DpsiIdT * dTdy + DpsiIdP * dPdy;
        	Dphi_im_dx[2] = DpsiIdR * dRdz + DpsiIdT * dTdz ;

                for( i=0; i<3; i++){
                *sourceSF1 += Dphi_re_dx[i] * Dphi_re_dx[i] + Dphi_im_dx[i] * Dphi_im_dx[i];
                }

                *sourceSF2 = mu * mu * ( psi_re * psi_re + psi_im * psi_im );
	
        }

	else {
		CCTK_WARN (0, "TwoPunctures_BBHSF: invalid choice for scalar field initial profile");
	}
}
