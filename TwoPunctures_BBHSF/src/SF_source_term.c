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
void SF_Initialize(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL *phi_re,
                   CCTK_REAL *phi_im) {
  /* Initialize the real and imaginary parts of the scalar field with an
   * analytic profile */

  DECLARE_CCTK_PARAMETERS;

  /*================ Define radius and angles ==============*/
  // Positions
  CCTK_REAL xp[3];
  xp[0] = x;
  xp[1] = y;
  xp[2] = z;

  // Coordinate radius rr and polar radial coordinate rho
  CCTK_REAL rr, rr2, irr, irr2;
  CCTK_REAL rho, rho2, irho, irho2;
  rr = sqrt(xp[0] * xp[0] + xp[1] * xp[1] + xp[2] * xp[2]);
  rho = sqrt(xp[0] * xp[0] + xp[1] * xp[1]);

  apply_cutoff(&rr); // Set lower bound for rr using TP_tiny
  apply_cutoff(&rho);
  rr2 = rr * rr;
  rho2 = rho * rho;

  irr = 1. / rr;
  irho = 1. / rho;
  irr2 = irr * irr;
  irho2 = irho * irho;

  CCTK_REAL deltR_w = (rr - r0SF) / widthSF;

  // Spherical polar coordinate angles
  CCTK_REAL ctheta, ctheta2, cphi, cphi2;
  CCTK_REAL stheta, stheta2, sphi, sphi2;

  ctheta = xp[2] * irr;
  stheta = rho * irr;
  cphi = xp[0] * irho;
  sphi = xp[1] * irho;

  ctheta2 = ctheta * ctheta;
  stheta2 = stheta * stheta;
  cphi2 = cphi * cphi;
  sphi2 = sphi * sphi;

  CCTK_REAL c2phi, s2phi;
  c2phi = cphi2 - sphi2;
  s2phi = 2.0 * cphi * sphi;

  /*========================================================*/
  // Currently implemented profiles:
  // Radially Gaussian with monopolar/dipolar angular dependence, centered at
  // the origin

  CCTK_INT monopole_statement, dipole_statement;

  /* Gaussian l=m=0 scalar field  */
  monopole_statement = (CCTK_EQUALS(initial_scalar, "ID_SF_Gaussian") &&
                        CCTK_EQUALS(scalar_GaussProfile, "single_mode") &&
                        (l0SF == 0) && (m0SF == 0));

  /* Gaussian l=m=1 scalar field  */
  dipole_statement = (CCTK_EQUALS(initial_scalar, "ID_SF_Gaussian") &&
                      CCTK_EQUALS(scalar_GaussProfile, "single_mode") &&
                      (l0SF == 1) && (m0SF == 1));

  if (monopole_statement) {
    *phi_re = (1. / sqrt(4 * Pi)) * ampSF * exp(-deltR_w * deltR_w);
    *phi_im = 0.;
  } else if (dipole_statement) {
    *phi_re =
        -ampSF * sqrt(1.5 / 4. / Pi) * stheta * cphi * exp(-deltR_w * deltR_w);
    *phi_im =
        -ampSF * sqrt(1.5 / 4. / Pi) * stheta * sphi * exp(-deltR_w * deltR_w);
  } else {
    CCTK_WARN(
        0,
        "TwoPunctures_BBHSF: invalid choice for scalar field initial profile");
  }
}

void conf_flat_analytic_SF_source_term(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
                                       CCTK_REAL phi_re, CCTK_REAL phi_im,
                                       CCTK_REAL *sourceSF1,
                                       CCTK_REAL *sourceSF2,
                                       CCTK_REAL Dphi_re_dx[3],
                                       CCTK_REAL Dphi_im_dx[3]) {
  /* Set the Cartesian derivatives and source terms from the scalar field
   * after the amplitude of the field was set using SF_initialize */

  DECLARE_CCTK_PARAMETERS;

  /*================ Define radius and angles ==============*/
  // Positions
  CCTK_REAL xp[3];
  xp[0] = x;
  xp[1] = y;
  xp[2] = z;

  // Coordinate radius rr and polar radial coordinate rho
  CCTK_REAL rr, rr2, irr, irr2;
  CCTK_REAL rho, rho2, irho, irho2;
  rr = sqrt(xp[0] * xp[0] + xp[1] * xp[1] + xp[2] * xp[2]);
  rho = sqrt(xp[0] * xp[0] + xp[1] * xp[1]);

  apply_cutoff(&rr); // Set lower bound for rr using TP_tiny
  apply_cutoff(&rho);
  rr2 = rr * rr;
  rho2 = rho * rho;

  irr = 1. / rr;
  irho = 1. / rho;
  irr2 = irr * irr;
  irho2 = irho * irho;

  CCTK_REAL deltR = (rr - r0SF);

  // Spherical polar coordinate angles
  CCTK_REAL ctheta, ctheta2, cphi, cphi2;
  CCTK_REAL stheta, stheta2, sphi, sphi2;

  ctheta = xp[2] * irr;
  stheta = rho * irr;
  cphi = xp[0] * irho;
  sphi = xp[1] * irho;

  ctheta2 = ctheta * ctheta;
  stheta2 = stheta * stheta;
  cphi2 = cphi * cphi;
  sphi2 = sphi * sphi;

  CCTK_REAL c2phi, s2phi;
  c2phi = cphi2 - sphi2;
  s2phi = 2.0 * cphi * sphi;

  /*--- Spherical coordinates to cartesian coordinates derivatives ---*/
  CCTK_REAL dRdx, dRdy, dRdz, dTdx, dTdy, dTdz, dPdx, dPdy;
  // dr/dx
  dRdx = xp[0] * irr;
  dRdy = xp[1] * irr;
  dRdz = xp[2] * irr;

  // d(theta)/dx
  dTdx = xp[0] * xp[2] * irho * irr2;
  dTdy = xp[1] * xp[2] * irho * irr2;
  dTdz = -rho * irr2;

  // d(phi)/dx
  dPdx = -xp[1] * irho2;
  dPdy = xp[0] * irho2;

  /*========================================================*/
  /* Analytically compute the derivatives of the scalar field components
   * Currently only supports scalar field profiles implemented in SF_Initialize
   */

  CCTK_INT monopole_statement, dipole_statement;

  /* Gaussian l=m=0 scalar field choice */
  monopole_statement = (CCTK_EQUALS(initial_scalar, "ID_SF_Gaussian") &&
                        CCTK_EQUALS(scalar_GaussProfile, "single_mode") &&
                        (l0SF == 0) && (m0SF == 0));

  /* Gaussian l=m=1 scalar field choice */
  dipole_statement = (CCTK_EQUALS(initial_scalar, "ID_SF_Gaussian") &&
                      CCTK_EQUALS(scalar_GaussProfile, "single_mode") &&
                      (l0SF == 1) && (m0SF == 1));

  /*--- Scalar field ---*/
  CCTK_REAL psi_re, psi_im;

  psi_re = phi_re;
  psi_im = phi_im;

  if (monopole_statement) {

    CCTK_REAL DpsiRdR;
    CCTK_REAL DpsiIdR;

    DpsiRdR = -2. * deltR * psi_re / (widthSF * widthSF);
    DpsiIdR = 0.;

    Dphi_re_dx[0] = DpsiRdR * dRdx;
    Dphi_re_dx[1] = DpsiRdR * dRdy;
    Dphi_re_dx[2] = DpsiRdR * dRdz;

    Dphi_im_dx[0] = DpsiIdR * dRdx;
    Dphi_im_dx[1] = DpsiIdR * dRdy;
    Dphi_im_dx[2] = DpsiIdR * dRdz;

  } else if (dipole_statement) {

    CCTK_REAL DpsiRdR, DpsiRdT, DpsiRdP;
    CCTK_REAL DpsiIdR, DpsiIdT, DpsiIdP;

    DpsiRdR = -2. * deltR * psi_re / (widthSF * widthSF);
    DpsiRdT = psi_re * ctheta / stheta;
    DpsiRdP = -psi_re * sphi / cphi;

    DpsiIdR = -2. * deltR * psi_im / (widthSF * widthSF);
    DpsiIdT = psi_im * ctheta / stheta;
    DpsiIdP = psi_im * cphi / sphi;

    Dphi_re_dx[0] = DpsiRdR * dRdx + DpsiRdT * dTdx + DpsiRdP * dPdx;
    Dphi_re_dx[1] = DpsiRdR * dRdy + DpsiRdT * dTdy + DpsiRdP * dPdy;
    Dphi_re_dx[2] = DpsiRdR * dRdz + DpsiRdT * dTdz;

    Dphi_im_dx[0] = DpsiIdR * dRdx + DpsiIdT * dTdx + DpsiIdP * dPdx;
    Dphi_im_dx[1] = DpsiIdR * dRdy + DpsiIdT * dTdy + DpsiIdP * dPdy;
    Dphi_im_dx[2] = DpsiIdR * dRdz + DpsiIdT * dTdz;

  } else {
    CCTK_WARN(
        0,
        "TwoPunctures_BBHSF: invalid choice for scalar field initial profile");
  }

  /*========================================================*/
  /* Compute source terms */

  // Kinetic term
  for (int i = 0; i < 3; i++) {
    *sourceSF1 += Dphi_re_dx[i] * Dphi_re_dx[i] + Dphi_im_dx[i] * Dphi_im_dx[i];
  }

  // Mass term
  *sourceSF2 = mu * mu * (psi_re * psi_re + psi_im * psi_im);
}
