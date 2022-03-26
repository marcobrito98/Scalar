/* TwoPuncturesSF:  File  "Equations.c"*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk_Parameters.h"
#include "TP_utilities.h"
#include "TwoPuncturesSF.h"

/* U.d0[ivar]   = U[ivar];  (ivar = 0..nvar-1) */
/* U.d1[ivar]   = U[ivar]_x;  */
/* U.d2[ivar]   = U[ivar]_y;  */
/* U.d3[ivar]   = U[ivar]_z;  */
/* U.d11[ivar]  = U[ivar]_xx; */
/* U.d12[ivar]  = U[ivar]_xy; */
/* U.d13[ivar]  = U[ivar]_xz;*/
/* U.d22[ivar]  = U[ivar]_yy;*/
/* U.d23[ivar]  = U[ivar]_yz;*/
/* U.d33[ivar]  = U[ivar]_zz;*/

CCTK_REAL
BY_KKofxyzSF (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z)
{
  DECLARE_CCTK_PARAMETERS;
  int i, j;
  CCTK_REAL r_plus, r2_plus, r3_plus, r_minus, r2_minus, r3_minus, np_Pp, nm_Pm,
    Aij, AijAij, n_plus[3], n_minus[3], np_Sp[3], nm_Sm[3];

  r2_plus = (x - par_b) * (x - par_b) + y * y + z * z;
  r2_minus = (x + par_b) * (x + par_b) + y * y + z * z;
  r_plus = sqrt (r2_plus);
  r_minus = sqrt (r2_minus);
  r3_plus = r_plus * r2_plus;
  r3_minus = r_minus * r2_minus;

  n_plus[0] = (x - par_b) / r_plus;
  n_minus[0] = (x + par_b) / r_minus;
  n_plus[1] = y / r_plus;
  n_minus[1] = y / r_minus;
  n_plus[2] = z / r_plus;
  n_minus[2] = z / r_minus;

  /* dot product: np_Pp = (n_+).(P_+); nm_Pm = (n_-).(P_-) */
  np_Pp = 0;
  nm_Pm = 0;
  for (i = 0; i < 3; i++)
  {
    np_Pp += n_plus[i] * par_P_plus[i];
    nm_Pm += n_minus[i] * par_P_minus[i];
  }
  /* cross product: np_Sp[i] = [(n_+) x (S_+)]_i; nm_Sm[i] = [(n_-) x (S_-)]_i*/
  np_Sp[0] = n_plus[1] * par_S_plus[2] - n_plus[2] * par_S_plus[1];
  np_Sp[1] = n_plus[2] * par_S_plus[0] - n_plus[0] * par_S_plus[2];
  np_Sp[2] = n_plus[0] * par_S_plus[1] - n_plus[1] * par_S_plus[0];
  nm_Sm[0] = n_minus[1] * par_S_minus[2] - n_minus[2] * par_S_minus[1];
  nm_Sm[1] = n_minus[2] * par_S_minus[0] - n_minus[0] * par_S_minus[2];
  nm_Sm[2] = n_minus[0] * par_S_minus[1] - n_minus[1] * par_S_minus[0];
  AijAij = 0;
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {				/* Bowen-York-Curvature :*/
      Aij =
	+ 1.5 * (par_P_plus[i] * n_plus[j] + par_P_plus[j] * n_plus[i]
                 + np_Pp * n_plus[i] * n_plus[j]) / r2_plus
	+ 1.5 * (par_P_minus[i] * n_minus[j] + par_P_minus[j] * n_minus[i]
		 + nm_Pm * n_minus[i] * n_minus[j]) / r2_minus
	- 3.0 * (np_Sp[i] * n_plus[j] + np_Sp[j] * n_plus[i]) / r3_plus
	- 3.0 * (nm_Sm[i] * n_minus[j] + nm_Sm[j] * n_minus[i]) / r3_minus;
      if (i == j)
	Aij -= +1.5 * (np_Pp / r2_plus + nm_Pm / r2_minus);
      AijAij += Aij * Aij;
    }
  }

  return AijAij;
}

void
BY_AijofxyzSF (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL Aij[3][3])
{
  DECLARE_CCTK_PARAMETERS;
  int i, j;
  CCTK_REAL r_plus, r2_plus, r3_plus, r_minus, r2_minus, r3_minus, np_Pp, nm_Pm,
    n_plus[3], n_minus[3], np_Sp[3], nm_Sm[3];

  r2_plus = (x - par_b) * (x - par_b) + y * y + z * z;
  r2_minus = (x + par_b) * (x + par_b) + y * y + z * z;
  r2_plus = sqrt (pow (r2_plus, 2) + pow (TP_epsilon, 4));
  r2_minus = sqrt (pow (r2_minus, 2) + pow (TP_epsilon, 4));
  if (r2_plus < pow(TP_Tiny,2))
    r2_plus = pow(TP_Tiny,2);
  if (r2_minus < pow(TP_Tiny,2))
    r2_minus = pow(TP_Tiny,2);
  r_plus = sqrt (r2_plus);
  r_minus = sqrt (r2_minus);
  r3_plus = r_plus * r2_plus;
  r3_minus = r_minus * r2_minus;

  n_plus[0] = (x - par_b) / r_plus;
  n_minus[0] = (x + par_b) / r_minus;
  n_plus[1] = y / r_plus;
  n_minus[1] = y / r_minus;
  n_plus[2] = z / r_plus;
  n_minus[2] = z / r_minus;

  /* dot product: np_Pp = (n_+).(P_+); nm_Pm = (n_-).(P_-) */
  np_Pp = 0;
  nm_Pm = 0;
  for (i = 0; i < 3; i++)
  {
    np_Pp += n_plus[i] * par_P_plus[i];
    nm_Pm += n_minus[i] * par_P_minus[i];
  }
  /* cross product: np_Sp[i] = [(n_+) x (S_+)]_i; nm_Sm[i] = [(n_-) x (S_-)]_i*/
  np_Sp[0] = n_plus[1] * par_S_plus[2] - n_plus[2] * par_S_plus[1];
  np_Sp[1] = n_plus[2] * par_S_plus[0] - n_plus[0] * par_S_plus[2];
  np_Sp[2] = n_plus[0] * par_S_plus[1] - n_plus[1] * par_S_plus[0];
  nm_Sm[0] = n_minus[1] * par_S_minus[2] - n_minus[2] * par_S_minus[1];
  nm_Sm[1] = n_minus[2] * par_S_minus[0] - n_minus[0] * par_S_minus[2];
  nm_Sm[2] = n_minus[0] * par_S_minus[1] - n_minus[1] * par_S_minus[0];
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {				/* Bowen-York-Curvature :*/
      Aij[i][j] =
        + 1.5 * (par_P_plus[i] * n_plus[j] + par_P_plus[j] * n_plus[i]
		 + np_Pp * n_plus[i] * n_plus[j]) / r2_plus
	+ 1.5 * (par_P_minus[i] * n_minus[j] + par_P_minus[j] * n_minus[i]
		 + nm_Pm * n_minus[i] * n_minus[j]) / r2_minus
	- 3.0 * (np_Sp[i] * n_plus[j] + np_Sp[j] * n_plus[i]) / r3_plus
	- 3.0 * (nm_Sm[i] * n_minus[j] + nm_Sm[j] * n_minus[i]) / r3_minus;
      if (i == j)
	Aij[i][j] -= +1.5 * (np_Pp / r2_plus + nm_Pm / r2_minus);
    }
  }
}

double dot(CCTK_REAL u[3], CCTK_REAL v[3]) {
    return  u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

/*-----------------------------------------------------------*/
/********           Nonlinear Equations                ***********/
/*-----------------------------------------------------------*/
void
NonLinEquationsSF (CCTK_REAL rho_adm,
     CCTK_REAL A, CCTK_REAL B, CCTK_REAL X, CCTK_REAL R,
		 CCTK_REAL x, CCTK_REAL r, CCTK_REAL phi,
		 CCTK_REAL y, CCTK_REAL z, derivs U, CCTK_REAL *values)
{
  DECLARE_CCTK_PARAMETERS;
  CCTK_REAL r_plus, r_minus, psi, psi2, psi4, psi7;
  CCTK_REAL sourceSF;

  r_plus    = sqrt( (x - par_b)*(x - par_b) + y*y + z*z );
  r_minus   = sqrt( (x + par_b)*(x + par_b) + y*y + z*z );

  psi   = 1. + 0.5 * par_m_plus / r_plus + 0.5 * par_m_minus / r_minus + U.d0[0];
  psi2  = psi * psi;
  psi4  = psi2 * psi2;
  psi7  = psi * psi2 * psi4;
  
  double dpsi[3];

  dpsi[0] =   - 0.5 * par_m_plus  / pow(r_plus , 3) * (x - par_b)
              - 0.5 * par_m_minus / pow(r_minus, 3) * (x + par_b)
              + U.d1[0];
  dpsi[1] =   - 0.5 * par_m_plus  / pow(r_plus , 3) * y
              - 0.5 * par_m_minus / pow(r_minus, 3) * y
              + U.d2[0];
  dpsi[2] =   - 0.5 * par_m_plus  / pow(r_plus , 3) * z
              - 0.5 * par_m_minus / pow(r_minus, 3) * z
              + U.d3[0];

  double rr;
  double Phibar1, dPhibar1[3];
  double Phibar2, dPhibar2[3];

  rr = sqrt( x*x + y*y + z*z);

  if (l0SF == 0) if (m0SF == 0) {
      Phibar1       = sqrt(0.25/Pi) * ampSF * exp( -pow((rr - r0SF)/widthSF, 2) ); 
      Phibar2       = 0; 

      dPhibar1[0]   = -2 * Phibar1 * (rr - r0SF) / pow(widthSF, 2) * x / rr;
      dPhibar1[1]   = -2 * Phibar1 * (rr - r0SF) / pow(widthSF, 2) * y / rr;
      dPhibar1[2]   = -2 * Phibar1 * (rr - r0SF) / pow(widthSF, 2) * z / rr;
      dPhibar2[0]   = -2 * Phibar2 * (rr - r0SF) / pow(widthSF, 2) * x / rr;
      dPhibar2[1]   = -2 * Phibar2 * (rr - r0SF) / pow(widthSF, 2) * y / rr;
      dPhibar2[2]   = -2 * Phibar2 * (rr - r0SF) / pow(widthSF, 2) * z / rr;
  }
  if (l0SF == 1) if (m0SF == 1) {
      double rho;
      double sthe, dsthe[3];
      double sphi, dsphi[3];
      double cphi, dcphi[3];

      rho       = sqrt( x*x + y*y );
      sthe      = rho / rr;
      sphi      = y / rho;
      cphi      = x / rho;

      dsthe[0]  = x / rho / rr - x * rho / pow(rr, 3);
      dsthe[1]  = y / rho / rr - y * rho / pow(rr, 3);
      dsthe[2]  = 0 / rho / rr - z * rho / pow(rr, 3);
      dsphi[0]  = 0 / rho - 1 * y / pow(rho, 3);
      dsphi[1]  = 1 / rho - 1 * y / pow(rho, 3);
      dsphi[2]  = 0 / rho - 0 * y / pow(rho, 3);
      dcphi[0]  = 1 / rho - 1 * x / pow(rho, 3);
      dcphi[1]  = 0 / rho - 1 * x / pow(rho, 3);
      dcphi[2]  = 0 / rho - 0 * x / pow(rho, 3);

      Phibar1   = -sqrt(0.375/Pi) * ampSF * exp( -pow((rr - r0SF)/widthSF, 2) ) * sthe * cphi; 
      Phibar2   = -sqrt(0.375/Pi) * ampSF * exp( -pow((rr - r0SF)/widthSF, 2) ) * sthe * sphi; 

      dPhibar1[0] = -2 * Phibar1 * (rr - r0SF) / pow(widthSF, 2) * x / rr
          + Phibar1 / sthe * dsthe[0]
          + Phibar1 / cphi * dcphi[0];
      dPhibar1[1] = -2 * Phibar1 * (rr - r0SF) / pow(widthSF, 2) * y / rr
          + Phibar1 / sthe * dsthe[1]
          + Phibar1 / cphi * dcphi[1];
      dPhibar1[2] = -2 * Phibar1 * (rr - r0SF) / pow(widthSF, 2) * z / rr
          + Phibar1 / sthe * dsthe[2]
          + Phibar1 / cphi * dcphi[2];

      dPhibar2[0] = -2 * Phibar2 * (rr - r0SF) / pow(widthSF, 2) * x / rr
          + Phibar2 / sthe * dsthe[0]
          + Phibar2 / sphi * dsphi[0];
      dPhibar2[1] = -2 * Phibar2 * (rr - r0SF) / pow(widthSF, 2) * y / rr
          + Phibar2 / sthe * dsthe[1]
          + Phibar2 / sphi * dsphi[1];
      dPhibar2[2] = -2 * Phibar2 * (rr - r0SF) / pow(widthSF, 2) * z / rr
          + Phibar2 / sthe * dsthe[2]
          + Phibar2 / sphi * dsphi[2];
  }

  values[0] =
    U.d11[0] + U.d22[0] + U.d33[0]
    + 0.125 * BY_KKofxyzSF (x, y, z) / psi7
    + Pi * pow(psi, 2*npsi+5) * mu * mu
        * ( Phibar1*Phibar1 + Phibar2*Phibar2 )
    + Pi * pow(psi, 2*npsi+1)
        * ( dot(dPhibar1, dPhibar1) + dot(dPhibar2, dPhibar2) )
    + Pi * pow(psi, 2*npsi+0) * npsi
        * (   Phibar1 * dot(dPhibar1, dpsi)
            + Phibar2 * dot(dPhibar2, dpsi) )
    + Pi * pow(psi, 2*npsi-1) * npsi * npsi
        * dot(dpsi, dpsi)
        * ( Phibar1*Phibar1 + Phibar2*Phibar2 )
    ;
}

/*-----------------------------------------------------------*/
/********               Linear Equations                ***********/
/*-----------------------------------------------------------*/
void
LinEquationsSF (CCTK_REAL A, CCTK_REAL B, CCTK_REAL X, CCTK_REAL R,
	      CCTK_REAL x, CCTK_REAL r, CCTK_REAL phi,
	      CCTK_REAL y, CCTK_REAL z, derivs dU, derivs U, CCTK_REAL *values)
{
  DECLARE_CCTK_PARAMETERS;
  CCTK_REAL r_plus, r_minus, psi, psi2, psi4, psi8;
  CCTK_REAL sourceSF;

  r_plus    = sqrt( (x - par_b)*(x - par_b) + y*y + z*z);
  r_minus   = sqrt( (x + par_b)*(x + par_b) + y*y + z*z);

  psi   = 1. + 0.5 * par_m_plus / r_plus + 0.5 * par_m_minus / r_minus + U.d0[0];
  psi2  = psi * psi;
  psi4  = psi2 * psi2;
  psi8  = psi4 * psi4;

  double dpsi[3];

  dpsi[0] =   - 0.5 * par_m_plus  / pow(r_plus , 3) * (x - par_b)
              - 0.5 * par_m_minus / pow(r_minus, 3) * (x + par_b)
              + U.d1[0];
  dpsi[1] =   - 0.5 * par_m_plus  / pow(r_plus , 3) * y
              - 0.5 * par_m_minus / pow(r_minus, 3) * y
              + U.d2[0];
  dpsi[2] =   - 0.5 * par_m_plus  / pow(r_plus , 3) * z
              - 0.5 * par_m_minus / pow(r_minus, 3) * z
              + U.d3[0];

  double rr;
  double Phibar1, dPhibar1[3];
  double Phibar2, dPhibar2[3];

  rr = sqrt( x*x + y*y + z*z);

  if (l0SF == 0) if (m0SF == 0) {
      Phibar1       = sqrt(0.25/Pi) * ampSF * exp( -pow((rr - r0SF)/widthSF, 2) ); 
      Phibar2       = 0; 

      dPhibar1[0]   = -2 * Phibar1 * (rr - r0SF) / pow(widthSF, 2) * x / rr;
      dPhibar1[1]   = -2 * Phibar1 * (rr - r0SF) / pow(widthSF, 2) * y / rr;
      dPhibar1[2]   = -2 * Phibar1 * (rr - r0SF) / pow(widthSF, 2) * z / rr;
      dPhibar2[0]   = -2 * Phibar2 * (rr - r0SF) / pow(widthSF, 2) * x / rr;
      dPhibar2[1]   = -2 * Phibar2 * (rr - r0SF) / pow(widthSF, 2) * y / rr;
      dPhibar2[2]   = -2 * Phibar2 * (rr - r0SF) / pow(widthSF, 2) * z / rr;
  }

  if (l0SF == 1) if (m0SF == 1) {
      double rho;
      double sthe, dsthe[3];
      double sphi, dsphi[3];
      double cphi, dcphi[3];

      rho = sqrt( x*x + y*y );
      if (rho == 0) {
          sthe  = 0;
          sphi  = 0;
          cphi  = 0;
      } else {
          sthe  = rho / rr;
          sphi  = y / rho;
          cphi  = x / rho;
      }

      dsthe[0]  = x / rho / rr - x * rho / pow(rr, 3);
      dsthe[1]  = y / rho / rr - y * rho / pow(rr, 3);
      dsthe[2]  = 0 / rho / rr - z * rho / pow(rr, 3);
      dsphi[0]  = 0 / rho - 1 * y / pow(rho, 3);
      dsphi[1]  = 1 / rho - 1 * y / pow(rho, 3);
      dsphi[2]  = 0 / rho - 0 * y / pow(rho, 3);
      dcphi[0]  = 1 / rho - 1 * x / pow(rho, 3);
      dcphi[1]  = 0 / rho - 1 * x / pow(rho, 3);
      dcphi[2]  = 0 / rho - 0 * x / pow(rho, 3);

      Phibar1   = -sqrt(0.375/Pi) * ampSF * exp( -pow((rr - r0SF)/widthSF, 2) ) * sthe * cphi; 
      Phibar2   = -sqrt(0.375/Pi) * ampSF * exp( -pow((rr - r0SF)/widthSF, 2) ) * sthe * sphi; 

      dPhibar1[0] = -2 * Phibar1 * (rr - r0SF) / pow(widthSF, 2) * x / rr
          + Phibar1 / sthe * dsthe[0]
          + Phibar1 / cphi * dcphi[0];
      dPhibar1[1] = -2 * Phibar1 * (rr - r0SF) / pow(widthSF, 2) * y / rr
          + Phibar1 / sthe * dsthe[1]
          + Phibar1 / cphi * dcphi[1];
      dPhibar1[2] = -2 * Phibar1 * (rr - r0SF) / pow(widthSF, 2) * z / rr
          + Phibar1 / sthe * dsthe[2]
          + Phibar1 / cphi * dcphi[2];

      dPhibar2[0] = -2 * Phibar2 * (rr - r0SF) / pow(widthSF, 2) * x / rr
          + Phibar2 / sthe * dsthe[0]
          + Phibar2 / sphi * dsphi[0];
      dPhibar2[1] = -2 * Phibar2 * (rr - r0SF) / pow(widthSF, 2) * y / rr
          + Phibar2 / sthe * dsthe[1]
          + Phibar2 / sphi * dsphi[1];
      dPhibar2[2] = -2 * Phibar2 * (rr - r0SF) / pow(widthSF, 2) * z / rr
          + Phibar2 / sthe * dsthe[2]
          + Phibar2 / sphi * dsphi[2];
  }

  values[0] = dU.d11[0] + dU.d22[0] + dU.d33[0]
    - 0.875 * BY_KKofxyzSF (x, y, z) / psi8 * dU.d0[0]
    + Pi * pow(psi, 2*npsi+4) * (2*npsi+5) * dU.d0[0] * mu * mu
        * ( Phibar1*Phibar1 + Phibar2*Phibar2 )
    + Pi * pow(psi, 2*npsi+0) * (2*npsi+1) * dU.d0[0]
        * ( dot(dPhibar1, dPhibar1) + dot(dPhibar2, dPhibar2) )
    + Pi * pow(psi, 2*npsi-1) * (2*npsi+0) * dU.d0[0] * npsi
        * (   Phibar1 * dot(dPhibar1, dpsi)
            + Phibar2 * dot(dPhibar2, dpsi) )
    + Pi * pow(psi, 2*npsi-2) * (2*npsi-1) * dU.d0[0] * npsi * npsi
        * ( Phibar1*Phibar1 + Phibar2*Phibar2 )
        * dot(dpsi, dpsi)

    + Pi * pow(psi, 2*npsi+0) * 2 * npsi
        * (   dU.d1[0] * ( Phibar1 * dPhibar1[0] + Phibar2 * dPhibar2[0] )
            + dU.d2[0] * ( Phibar1 * dPhibar1[1] + Phibar2 * dPhibar2[1] )
            + dU.d3[0] * ( Phibar1 * dPhibar1[2] + Phibar2 * dPhibar2[2] )  )
    + Pi * pow(psi, 2*npsi-1) * 2 * npsi * npsi
        * ( Phibar1*Phibar1 + Phibar2*Phibar2 )
        * (   dU.d1[0] * dpsi[0]
            + dU.d2[0] * dpsi[1]
            + dU.d3[0] * dpsi[2] )
    ;
}

/*-----------------------------------------------------------*/
