/* TwoPunctures:  File  "Equations.c"*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk_Parameters.h"
#include "TP_utilities.h"
#include "TwoPunctures_BBHSF.h"

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
BBHSF_BY_KKofxyz(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z) {
  DECLARE_CCTK_PARAMETERS;
  int i, j;
  CCTK_REAL r_plus, r2_plus, r3_plus, r_minus, r2_minus, r3_minus, np_Pp, nm_Pm,
      Aij, AijAij, n_plus[3], n_minus[3], np_Sp[3], nm_Sm[3];

  r2_plus = (x - par_b) * (x - par_b) + y * y + z * z;
  r2_minus = (x + par_b) * (x + par_b) + y * y + z * z;
  r_plus = sqrt(r2_plus);
  r_minus = sqrt(r2_minus);
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
  for (i = 0; i < 3; i++) {
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
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) { /* Bowen-York-Curvature :*/
      Aij = +1.5 *
                (par_P_plus[i] * n_plus[j] + par_P_plus[j] * n_plus[i] +
                 np_Pp * n_plus[i] * n_plus[j]) /
                r2_plus +
            1.5 *
                (par_P_minus[i] * n_minus[j] + par_P_minus[j] * n_minus[i] +
                 nm_Pm * n_minus[i] * n_minus[j]) /
                r2_minus -
            3.0 * (np_Sp[i] * n_plus[j] + np_Sp[j] * n_plus[i]) / r3_plus -
            3.0 * (nm_Sm[i] * n_minus[j] + nm_Sm[j] * n_minus[i]) / r3_minus;

      if (i == j)
        Aij -= +1.5 * (np_Pp / r2_plus + nm_Pm / r2_minus);

      AijAij += Aij * Aij;
    }
  }

  return AijAij;
}

void BBHSF_BY_Aijofxyz(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
                       CCTK_REAL Aij[3][3]) {
  DECLARE_CCTK_PARAMETERS;
  int i, j;
  CCTK_REAL r_plus, r2_plus, r3_plus, r_minus, r2_minus, r3_minus, np_Pp, nm_Pm,
      n_plus[3], n_minus[3], np_Sp[3], nm_Sm[3];

  r2_plus = (x - par_b) * (x - par_b) + y * y + z * z;
  r2_minus = (x + par_b) * (x + par_b) + y * y + z * z;
  r2_plus = sqrt(pow(r2_plus, 2) + pow(TP_epsilon, 4));
  r2_minus = sqrt(pow(r2_minus, 2) + pow(TP_epsilon, 4));
  if (r2_plus < pow(TP_Tiny, 2))
    r2_plus = pow(TP_Tiny, 2);
  if (r2_minus < pow(TP_Tiny, 2))
    r2_minus = pow(TP_Tiny, 2);

  r_plus = sqrt(r2_plus);
  r_minus = sqrt(r2_minus);
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
  for (i = 0; i < 3; i++) {
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
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) { /* Bowen-York-Curvature :*/
      Aij[i][j] =
          +1.5 *
              (par_P_plus[i] * n_plus[j] + par_P_plus[j] * n_plus[i] +
               np_Pp * n_plus[i] * n_plus[j]) /
              r2_plus +
          1.5 *
              (par_P_minus[i] * n_minus[j] + par_P_minus[j] * n_minus[i] +
               nm_Pm * n_minus[i] * n_minus[j]) /
              r2_minus -
          3.0 * (np_Sp[i] * n_plus[j] + np_Sp[j] * n_plus[i]) / r3_plus -
          3.0 * (nm_Sm[i] * n_minus[j] + nm_Sm[j] * n_minus[i]) / r3_minus;
      if (i == j)
        Aij[i][j] -= +1.5 * (np_Pp / r2_plus + nm_Pm / r2_minus);
    }
  }
}

/*-----------------------------------------------------------*/
/********           Nonlinear Equations                ***********/
/*-----------------------------------------------------------*/
void BBHSF_NonLinEquations(CCTK_REAL rho_adm, CCTK_REAL A, CCTK_REAL B,
                           CCTK_REAL X, CCTK_REAL R, CCTK_REAL x, CCTK_REAL r,
                           CCTK_REAL phi, CCTK_REAL y, CCTK_REAL z, derivs U,
                           CCTK_REAL *values) {
  DECLARE_CCTK_PARAMETERS;
  // rho_adm is not supported

  // Distance to plus and minus punctures
  CCTK_REAL r_plus, r_minus;
  r_plus = sqrt((x - par_b) * (x - par_b) + y * y + z * z);
  r_minus = sqrt((x + par_b) * (x + par_b) + y * y + z * z);
  apply_cutoff(&r_plus);
  apply_cutoff(&r_minus);

  // Inverse of distance to plus and minus punctures
  CCTK_REAL ir_plus, ir_minus;
  ir_plus = 1. / r_plus;
  ir_minus = 1. / r_minus;

  CCTK_REAL irplus3, irminus3;
  irplus3 = ir_plus * ir_plus * ir_plus;
  irminus3 = ir_minus * ir_minus * ir_minus;

  /*========================================================*/
  /* Conformal factor */

  // Write conformal factor as Brill-Lindquist solution (psiBL) + correction u
  CCTK_REAL psi, psi2, psi4, psi7;
  psi =
      1. + 0.5 * par_m_plus * ir_plus + 0.5 * par_m_minus * ir_minus + U.d0[0];
  psi2 = psi * psi;
  psi4 = psi2 * psi2;
  psi7 = psi4 * psi2 * psi;

  // Derivatives of psiBL
  CCTK_REAL DpsiBLdx[3];
  DpsiBLdx[0] = -0.5 * par_m_plus * (x - par_b) * irplus3 -
                0.5 * par_m_minus * (x + par_b) * irminus3;
  DpsiBLdx[1] =
      -0.5 * par_m_plus * y * irplus3 - 0.5 * par_m_minus * y * irminus3;
  DpsiBLdx[2] =
      -0.5 * par_m_plus * z * irplus3 - 0.5 * par_m_minus * z * irminus3;

  // Derivatives of psi
  CCTK_REAL Dpsidx[3];
  Dpsidx[0] = DpsiBLdx[0] + U.d1[0];
  Dpsidx[1] = DpsiBLdx[1] + U.d2[0];
  Dpsidx[2] = DpsiBLdx[2] + U.d3[0];

  /*========================================================*/
  /* Scalar field */

  // Initialize the real and imaginary parts of the scalar field and its
  // derivatives
  CCTK_REAL phi_re, phi_im;
  CCTK_REAL Dphi_re_dx[3], Dphi_im_dx[3];

  phi_re = 0.;
  phi_im = 0.;
  Dphi_re_dx[0] = Dphi_re_dx[1] = Dphi_re_dx[2] = 0.;
  Dphi_im_dx[0] = Dphi_im_dx[1] = Dphi_im_dx[2] = 0.;

  // Contributions of the scalar field to the Hamiltonian constraint
  CCTK_REAL sourceSF1, sourceSF2;
  sourceSF1 = 0.;
  sourceSF2 = 0.;

  // With backreaction on, we use these functions to fill in the scalar field
  // itself and its contributions to H
  if (switch_on_backreaction) {

    // Set real and imag parts of phi
    SF_Initialize(x, y, z, &phi_re, &phi_im);

    // Compute Cartesian derivatives of phi and source terms
    conf_flat_analytic_SF_source_term(x, y, z, phi_re, phi_im, &sourceSF1,
                                      &sourceSF2, &Dphi_re_dx, &Dphi_im_dx);
  }

  /*========================================================*/
  /* Compute the terms in the linearized Hamiltonian constraint */
  CCTK_REAL term1, term2, term3, term4, term5, term6;

  // Flat laplacian of u
  term1 = U.d11[0] + U.d22[0] + U.d33[0];

  // Bowen-York curvature
  term2 = 0.125 * BBHSF_BY_KKofxyz(x, y, z) / psi7;

  // Scalar field spatial derivatives
  term3 = Pi * sourceSF1 * pow(psi, 2 * delta + 1);

  // Scalar field mass term
  term4 = Pi * sourceSF2 * pow(psi, 2 * delta + 5);

  // Mixed psi and field derivatives term
  term5 = 2 * Pi * delta * pow(psi, 2 * delta) *
          (phi_re * scalarproduct(Dpsidx, Dphi_re_dx, 3) +
           phi_im * scalarproduct(Dpsidx, Dphi_im_dx, 3));

  // Psi derivatives term
  term6 = Pi * delta * delta * pow(psi, 2 * delta - 1) *
          (phi_re * phi_re + phi_im * phi_im) *
          scalarproduct(Dpsidx, Dpsidx, 3);

  // Hamiltonian constraint equation
  values[0] = term1 + term2 + term3 + term4 + term5 + term6;
}

/*-----------------------------------------------------------*/
/********               Linear Equations                ***********/
/*-----------------------------------------------------------*/
void BBHSF_LinEquations(CCTK_REAL A, CCTK_REAL B, CCTK_REAL X, CCTK_REAL R,
                        CCTK_REAL x, CCTK_REAL r, CCTK_REAL phi, CCTK_REAL y,
                        CCTK_REAL z, derivs dU, derivs U, CCTK_REAL *values) {
  DECLARE_CCTK_PARAMETERS;

  // Distance to plus and minus punctures
  CCTK_REAL r_plus, r_minus;
  r_plus = sqrt((x - par_b) * (x - par_b) + y * y + z * z);
  r_minus = sqrt((x + par_b) * (x + par_b) + y * y + z * z);
  apply_cutoff(&r_plus);
  apply_cutoff(&r_minus);

  // Inverse of distance to plus and minus punctures
  CCTK_REAL ir_plus, ir_minus;
  ir_plus = 1. / r_plus;
  ir_minus = 1. / r_minus;

  CCTK_REAL irplus3, irminus3;
  irplus3 = ir_plus * ir_plus * ir_plus;
  irminus3 = ir_minus * ir_minus * ir_minus;

  /*========================================================*/
  /* Conformal factor */

  // Write conformal factor as Brill-Lindquist solution (psiBL) + correction u
  CCTK_REAL psi, psi2, psi4, psi8;
  psi =
      1. + 0.5 * par_m_plus * ir_plus + 0.5 * par_m_minus * ir_minus + U.d0[0];
  psi2 = psi * psi;
  psi4 = psi2 * psi2;
  psi8 = psi4 * psi4;

  // Derivatives of psiBL
  CCTK_REAL DpsiBLdx[3];
  DpsiBLdx[0] = -0.5 * par_m_plus * (x - par_b) * irplus3 -
                0.5 * par_m_minus * (x + par_b) * irminus3;
  DpsiBLdx[1] =
      -0.5 * par_m_plus * y * irplus3 - 0.5 * par_m_minus * y * irminus3;
  DpsiBLdx[2] =
      -0.5 * par_m_plus * z * irplus3 - 0.5 * par_m_minus * z * irminus3;

  // Derivatives of psi
  CCTK_REAL Dpsidx[3];
  Dpsidx[0] = DpsiBLdx[0] + U.d1[0];
  Dpsidx[1] = DpsiBLdx[1] + U.d2[0];
  Dpsidx[2] = DpsiBLdx[2] + U.d3[0];

  /*========================================================*/
  /* Scalar field */

  // Initialize the real and imaginary parts of the scalar field and its
  // derivatives
  CCTK_REAL phi_re, phi_im;
  CCTK_REAL Dphi_re_dx[3], Dphi_im_dx[3];

  phi_re = 0.;
  phi_im = 0.;
  Dphi_re_dx[0] = Dphi_re_dx[1] = Dphi_re_dx[2] = 0.;
  Dphi_im_dx[0] = Dphi_im_dx[1] = Dphi_im_dx[2] = 0.;

  // Contributions of the scalar field to the Hamiltonian constraint
  CCTK_REAL sourceSF1, sourceSF2;
  sourceSF1 = 0.;
  sourceSF2 = 0.;

  // With backreaction on, we use these functions to fill in the scalar field
  // itself and its contributions to H
  if (switch_on_backreaction) {

    // Set real and imag parts of phi
    SF_Initialize(x, y, z, &phi_re, &phi_im);

    // Compute Cartesian derivatives of phi and source terms
    conf_flat_analytic_SF_source_term(x, y, z, phi_re, phi_im, &sourceSF1,
                                      &sourceSF2, &Dphi_re_dx, &Dphi_im_dx);
  }

  /*========================================================*/
  /* Compute the terms in the linearized Hamiltonian constraint */

  CCTK_REAL term1, term2, term3, term4, term5, term6;

  // Flat laplacian of linearized u
  term1 = dU.d11[0] + dU.d22[0] + dU.d33[0];

  // Bowen-York curvature linearized term
  term2 = -0.875 * BBHSF_BY_KKofxyz(x, y, z) / psi8 * dU.d0[0];

  // Scalar field spatial derivatives linearized term
  term3 = Pi * (2 * delta + 1) * pow(psi, 2 * delta) * sourceSF1 * dU.d0[0];

  // Scalar field mass linearized term
  term4 = Pi * (2 * delta + 5) * pow(psi, 2 * delta + 4) * sourceSF2 * dU.d0[0];

  // Mixed psi and field derivatives linearized term
  term5 = 2 * Pi * delta *
          (phi_re * (pow(psi, 2 * delta) *
                         (dU.d1[0] * Dphi_re_dx[0] + dU.d2[0] * Dphi_re_dx[1] +
                          dU.d3[0] * Dphi_re_dx[2]) +
                     pow(psi, 2 * delta - 1) * 2 * delta * dU.d0[0] *
                         scalarproduct(Dpsidx, Dphi_re_dx, 3)) +
           phi_im * (pow(psi, 2 * delta) *
                         (dU.d1[0] * Dphi_im_dx[0] + dU.d2[0] * Dphi_im_dx[1] +
                          dU.d3[0] * Dphi_im_dx[2]) +
                     pow(psi, 2 * delta - 1) * 2 * delta * dU.d0[0] *
                         scalarproduct(Dpsidx, Dphi_im_dx, 3)));

  // Psi derivatives linearized term
  term6 = Pi * delta * delta * (phi_re * phi_re + phi_im * phi_im) *
          (pow(psi, 2 * delta - 1) * 2 *
               (dU.d1[0] * Dpsidx[0] + dU.d2[0] * Dpsidx[1] +
                dU.d3[0] * Dpsidx[2]) +
           pow(psi, 2 * delta - 2) * (2 * delta - 1) * dU.d0[0] *
               scalarproduct(Dpsidx, Dpsidx, 3));

  // Linearized Hamiltonian constraint equation
  values[0] = term1 + term2 + term3 + term4 + term5 + term6;
}

/*-----------------------------------------------------------*/
