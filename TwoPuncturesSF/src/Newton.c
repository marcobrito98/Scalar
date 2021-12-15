/* TwoPuncturesSF:  File  "Newton.c"*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include "cctk_Parameters.h"
#include "TP_utilities.h"
#include "TwoPuncturesSF.h"

static int
bicgstabSF (CCTK_POINTER_TO_CONST const cctkGH,
          int const nvar, int const n1, int const n2, int const n3,
          derivs v, derivs dv,
          int const output, int const itmax, CCTK_REAL const tol,
          CCTK_REAL * restrict const normres);
static CCTK_REAL
norm_infSF (CCTK_REAL const * restrict const F,
          int const ntotal);
static void
relaxSF (CCTK_REAL * restrict const dv,
       int const nvar, int const n1, int const n2, int const n3,
       CCTK_REAL const * restrict const rhs,
       int const * restrict const ncols,
       int const * restrict const * restrict const cols,
       CCTK_REAL const * restrict const * restrict const JFD);
static void
residSF (CCTK_REAL * restrict const res,
       int const ntotal,
       CCTK_REAL const * restrict const dv,
       CCTK_REAL const * restrict const rhs,
       int const * restrict const ncols,
       int const * restrict const * restrict const cols,
       CCTK_REAL const * restrict const * restrict const JFD);
static void
LineRelax_alSF (CCTK_REAL * restrict const dv,
              int const j, int const k, int const nvar,
              int const n1, int const n2, int const n3,
	      CCTK_REAL const * restrict const rhs,
              int const * restrict const ncols,
              int const * restrict const * restrict const cols,
              CCTK_REAL const * restrict const * restrict const JFD);
static void
LineRelax_beSF (CCTK_REAL * restrict const dv,
              int const i, int const k, int const nvar,
              int const n1, int const n2, int const n3,
	      CCTK_REAL const * restrict const rhs,
              int const * restrict const ncols,
              int const * restrict const * restrict const cols,
              CCTK_REAL const * restrict const * restrict const JFD);
/* --------------------------------------------------------------------------*/
static CCTK_REAL
norm_infSF (CCTK_REAL const * restrict const F,
          int const ntotal)
{
  CCTK_REAL dmax = -1;
#pragma omp parallel
  {
    CCTK_REAL dmax1 = -1;
#pragma omp for
    for (int j = 0; j < ntotal; j++)
      if (fabs (F[j]) > dmax1)
        dmax1 = fabs (F[j]);
#pragma omp critical
    if (dmax1 > dmax)
      dmax = dmax1;
  }
  return dmax;
}
/* --------------------------------------------------------------------------*/
static void
residSF (CCTK_REAL * restrict const res,
       int const ntotal,
       CCTK_REAL const * restrict const dv,
       CCTK_REAL const * restrict const rhs,
       int const * restrict const ncols,
       int const * restrict const * restrict const cols,
       CCTK_REAL const * restrict const * restrict const JFD)
{
#pragma omp parallel for
  for (int i = 0; i < ntotal; i++)
  {
    CCTK_REAL JFDdv_i = 0;
    for (int m = 0; m < ncols[i]; m++)
      JFDdv_i += JFD[i][m] * dv[cols[i][m]];
    res[i] = rhs[i] - JFDdv_i;
  }
}

/* -------------------------------------------------------------------------*/
static void
LineRelax_alSF (CCTK_REAL * restrict const dv,
              int const j, int const k, int const nvar,
              int const n1, int const n2, int const n3,
	      CCTK_REAL const * restrict const rhs,
              int const * restrict const ncols,
              int const * restrict const * restrict const cols,
              CCTK_REAL const * restrict const * restrict const JFD)
{
  int i, m, Ic, Ip, Im, col, ivar;

  gsl_vector *diag = gsl_vector_alloc(n1);
  gsl_vector *e = gsl_vector_alloc(n1-1); /* above diagonal */
  gsl_vector *f = gsl_vector_alloc(n1-1); /* below diagonal */
  gsl_vector *b = gsl_vector_alloc(n1);   /* rhs */
  gsl_vector *x = gsl_vector_alloc(n1);   /* solution vector */

  for (ivar = 0; ivar < nvar; ivar++)
  {
    gsl_vector_set_zero(diag);
    gsl_vector_set_zero(e);
    gsl_vector_set_zero(f);
    for (i = 0; i < n1; i++)
    {
      Ip = IndexSF (ivar, i + 1, j, k, nvar, n1, n2, n3);
      Ic = IndexSF (ivar, i, j, k, nvar, n1, n2, n3);
      Im = IndexSF (ivar, i - 1, j, k, nvar, n1, n2, n3);
      gsl_vector_set(b,i,rhs[Ic]);
      for (m = 0; m < ncols[Ic]; m++)
      {
	col = cols[Ic][m];
	if (col != Ip && col != Ic && col != Im)
          *gsl_vector_ptr(b, i) -= JFD[Ic][m] * dv[col];
	else
	{
	  if (col == Im && i > 0)
            gsl_vector_set(f,i-1,JFD[Ic][m]);
	  if (col == Ic)
            gsl_vector_set(diag,i,JFD[Ic][m]);
	  if (col == Ip && i < n1-1)
            gsl_vector_set(e,i,JFD[Ic][m]);
	}
      }
    }
    gsl_linalg_solve_tridiag(diag, e, f, b, x);
    for (i = 0; i < n1; i++)
    {
      Ic = IndexSF (ivar, i, j, k, nvar, n1, n2, n3);
      dv[Ic] = gsl_vector_get(x, i);
    }
  }

  gsl_vector_free(diag);
  gsl_vector_free(e);
  gsl_vector_free(f);
  gsl_vector_free(b);
  gsl_vector_free(x);
}

/* --------------------------------------------------------------------------*/
static void
LineRelax_beSF (CCTK_REAL * restrict const dv,
              int const i, int const k, int const nvar,
              int const n1, int const n2, int const n3,
	      CCTK_REAL const * restrict const rhs,
              int const * restrict const ncols,
              int const * restrict const * restrict const cols,
              CCTK_REAL const * restrict const * restrict const JFD)
{
  int j, m, Ic, Ip, Im, col, ivar;

  gsl_vector *diag = gsl_vector_alloc(n2);
  gsl_vector *e = gsl_vector_alloc(n2-1); /* above diagonal */
  gsl_vector *f = gsl_vector_alloc(n2-1); /* below diagonal */
  gsl_vector *b = gsl_vector_alloc(n2);   /* rhs */
  gsl_vector *x = gsl_vector_alloc(n2);   /* solution vector */

  for (ivar = 0; ivar < nvar; ivar++)
  {
    gsl_vector_set_zero(diag);
    gsl_vector_set_zero(e);
    gsl_vector_set_zero(f);
    for (j = 0; j < n2; j++)
    {
      Ip = IndexSF (ivar, i, j + 1, k, nvar, n1, n2, n3);
      Ic = IndexSF (ivar, i, j, k, nvar, n1, n2, n3);
      Im = IndexSF (ivar, i, j - 1, k, nvar, n1, n2, n3);
      gsl_vector_set(b,j,rhs[Ic]);
      for (m = 0; m < ncols[Ic]; m++)
      {
	col = cols[Ic][m];
	if (col != Ip && col != Ic && col != Im)
          *gsl_vector_ptr(b, j) -= JFD[Ic][m] * dv[col];
	else
	{
	  if (col == Im && j > 0)
            gsl_vector_set(f,j-1,JFD[Ic][m]);
	  if (col == Ic)
            gsl_vector_set(diag,j,JFD[Ic][m]);
	  if (col == Ip && j < n2-1)
            gsl_vector_set(e,j,JFD[Ic][m]);
	}
      }
    }
    gsl_linalg_solve_tridiag(diag, e, f, b, x);
    for (j = 0; j < n2; j++)
    {
      Ic = IndexSF (ivar, i, j, k, nvar, n1, n2, n3);
      dv[Ic] = gsl_vector_get(x, j);
    }
  }
  gsl_vector_free(diag);
  gsl_vector_free(e);
  gsl_vector_free(f);
  gsl_vector_free(b);
  gsl_vector_free(x);
}

/* --------------------------------------------------------------------------*/
static void
relaxSF (CCTK_REAL * restrict const dv,
       int const nvar, int const n1, int const n2, int const n3,
       CCTK_REAL const * restrict const rhs,
       int const * restrict const ncols,
       int const * restrict const * restrict const cols,
       CCTK_REAL const * restrict const * restrict const JFD)
{
  int i, j, k, n;

  for (k = 0; k < n3; k = k + 2)
  {
    for (n = 0; n < N_PlaneRelax; n++)
    {
#pragma omp parallel for schedule(dynamic)
      for (i = 2; i < n1; i = i + 2)
	LineRelax_beSF (dv, i, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
#pragma omp parallel for schedule(dynamic)
      for (i = 1; i < n1; i = i + 2)
	LineRelax_beSF (dv, i, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
#pragma omp parallel for schedule(dynamic)
      for (j = 1; j < n2; j = j + 2)
	LineRelax_alSF (dv, j, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
#pragma omp parallel for schedule(dynamic)
      for (j = 0; j < n2; j = j + 2)
	LineRelax_alSF (dv, j, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
    }
  }
  for (k = 1; k < n3; k = k + 2)
  {
    for (n = 0; n < N_PlaneRelax; n++)
    {
#pragma omp parallel for schedule(dynamic)
      for (i = 0; i < n1; i = i + 2)
	LineRelax_beSF (dv, i, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
#pragma omp parallel for schedule(dynamic)
      for (i = 1; i < n1; i = i + 2)
	LineRelax_beSF (dv, i, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
#pragma omp parallel for schedule(dynamic)
      for (j = 1; j < n2; j = j + 2)
	LineRelax_alSF (dv, j, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
#pragma omp parallel for schedule(dynamic)
      for (j = 0; j < n2; j = j + 2)
	LineRelax_alSF (dv, j, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
    }
  }
}

/* --------------------------------------------------------------------------*/
void
TestRelaxSF (CCTK_POINTER_TO_CONST cctkGH,
     int nvar, int n1, int n2, int n3, derivs v,
	   CCTK_REAL *dv)
{
  DECLARE_CCTK_PARAMETERS;
  int ntotal = n1 * n2 * n3 * nvar, **cols, *ncols, maxcol =
    StencilSize * nvar, j;
  CCTK_REAL *F, *res, **JFD;
  derivs u;

  F = dvectorSF (0, ntotal - 1);
  res = dvectorSF (0, ntotal - 1);
  allocate_derivs (&u, ntotal);

  JFD = dmatrixSF (0, ntotal - 1, 0, maxcol - 1);
  cols = imatrixSF (0, ntotal - 1, 0, maxcol - 1);
  ncols = ivectorSF (0, ntotal - 1);

  F_of_vSF (cctkGH, nvar, n1, n2, n3, v, F, u);

  SetMatrix_JFDSF (nvar, n1, n2, n3, u, ncols, cols, JFD);

  for (j = 0; j < ntotal; j++)
    dv[j] = 0;
  residSF (res, ntotal, dv, F, ncols, cols, JFD);
  printf ("Before: |F|=%20.15e\n", (double) norm1SF (res, ntotal));
  fflush(stdout);
  for (j = 0; j < NRELAX; j++)
  {
    relaxSF (dv, nvar, n1, n2, n3, F, ncols, cols, JFD);	/* solves JFD*sh = s*/
    if (j % Step_Relax == 0)
    {
      residSF (res, ntotal, dv, F, ncols, cols, JFD);
      printf ("j=%d\t |F|=%20.15e\n", j, (double) norm1SF (res, ntotal));
      fflush(stdout);
    }
  }

  residSF (res, ntotal, dv, F, ncols, cols, JFD);
  printf ("After: |F|=%20.15e\n", (double) norm1SF (res, ntotal));
  fflush(stdout);

  free_dvectorSF (F, 0, ntotal - 1);
  free_dvectorSF (res, 0, ntotal - 1);
  free_derivs (&u, ntotal);

  free_dmatrixSF (JFD, 0, ntotal - 1, 0, maxcol - 1);
  free_imatrixSF (cols, 0, ntotal - 1, 0, maxcol - 1);
  free_ivectorSF (ncols, 0, ntotal - 1);
}

/* --------------------------------------------------------------------------*/
static int
bicgstabSF (CCTK_POINTER_TO_CONST const cctkGH,
          int const nvar, int const n1, int const n2, int const n3,
          derivs v, derivs dv,
          int const output, int const itmax, CCTK_REAL const tol,
          CCTK_REAL * restrict const normres)
{
  DECLARE_CCTK_PARAMETERS;
  int ntotal = n1 * n2 * n3 * nvar, ii;
  CCTK_REAL alpha = 0, beta = 0;
  CCTK_REAL rho = 0, rho1 = 1, rhotol = 1e-50;
  CCTK_REAL omega = 0, omegatol = 1e-50;
  CCTK_REAL *p, *rt, *s, *t, *r, *vv;
  CCTK_REAL **JFD;
  int **cols, *ncols, maxcol = StencilSize * nvar;
  CCTK_REAL *F;
  derivs u, ph, sh;

  F = dvectorSF (0, ntotal - 1);
  allocate_derivs (&u, ntotal);

  JFD = dmatrixSF (0, ntotal - 1, 0, maxcol - 1);
  cols = imatrixSF (0, ntotal - 1, 0, maxcol - 1);
  ncols = ivectorSF (0, ntotal - 1);

  F_of_vSF (cctkGH, nvar, n1, n2, n3, v, F, u);
  SetMatrix_JFD (nvar, n1, n2, n3, u, ncols, cols, JFD);

  /* temporary storage */
  r = dvectorSF (0, ntotal - 1);
  p = dvectorSF (0, ntotal - 1);
  allocate_derivs (&ph, ntotal);
/*      ph  = dvectorSF(0, ntotal-1);*/
  rt = dvectorSF (0, ntotal - 1);
  s = dvectorSF (0, ntotal - 1);
  allocate_derivs (&sh, ntotal);
/*      sh  = dvectorSF(0, ntotal-1);*/
  t = dvectorSF (0, ntotal - 1);
  vv = dvectorSF (0, ntotal - 1);

  /* check */
  if (output == 1) {
    printf ("bicgstab:  itmax %d, tol %e\n", itmax, (double)tol);
    fflush(stdout);
  }

  /* compute initial residual rt = r = F - J*dv */
  J_times_dvSF (nvar, n1, n2, n3, dv, r, u);
#pragma omp parallel for
  for (int j = 0; j < ntotal; j++)
    rt[j] = r[j] = F[j] - r[j];

  *normres = norm2SF (r, ntotal);
  if (output == 1) {
    printf ("bicgstab: %5d  %10.3e\n", 0, (double) *normres);
    fflush(stdout);
  }

  if (*normres <= tol)
    return 0;

  /* cgs iteration */
  for (ii = 0; ii < itmax; ii++)
  {
    rho = scalarproductSF (rt, r, ntotal);
    if (fabs (rho) < rhotol)
      break;

    /* compute direction vector p */
    if (ii == 0)
    {
#pragma omp parallel for
      for (int j = 0; j < ntotal; j++)
	p[j] = r[j];
    }
    else
    {
      beta = (rho / rho1) * (alpha / omega);
#pragma omp parallel for
      for (int j = 0; j < ntotal; j++)
	p[j] = r[j] + beta * (p[j] - omega * vv[j]);
    }

    /* compute direction adjusting vector ph and scalar alpha */
#pragma omp parallel for
    for (int j = 0; j < ntotal; j++)
      ph.d0[j] = 0;
    for (int j = 0; j < NRELAX; j++)	/* solves JFD*ph = p by relaxation*/
      relaxSF (ph.d0, nvar, n1, n2, n3, p, ncols, cols, JFD);

    J_times_dvSF (nvar, n1, n2, n3, ph, vv, u);	/* vv=J*ph*/
    alpha = rho / scalarproductSF (rt, vv, ntotal);
#pragma omp parallel for
    for (int j = 0; j < ntotal; j++)
      s[j] = r[j] - alpha * vv[j];

    /* early check of tolerance */
    *normres = norm2SF (s, ntotal);
    if (*normres <= tol)
    {
#pragma omp parallel for
      for (int j = 0; j < ntotal; j++)
	dv.d0[j] += alpha * ph.d0[j];
      if (output == 1) {
	printf ("bicgstab: %5d  %10.3e  %10.3e  %10.3e  %10.3e\n",
		ii + 1, (double) *normres, (double)alpha, (double)beta, (double)omega);
        fflush(stdout);
      }
      break;
    }

    /* compute stabilizer vector sh and scalar omega */
#pragma omp parallel for
    for (int j = 0; j < ntotal; j++)
      sh.d0[j] = 0;
    for (int j = 0; j < NRELAX; j++)	/* solves JFD*sh = s by relaxation*/
      relaxSF (sh.d0, nvar, n1, n2, n3, s, ncols, cols, JFD);

    J_times_dvSF (nvar, n1, n2, n3, sh, t, u);	/* t=J*sh*/
    omega = scalarproductSF (t, s, ntotal) / scalarproductSF (t, t, ntotal);

    /* compute new solution approximation */
#pragma omp parallel for
    for (int j = 0; j < ntotal; j++)
    {
      dv.d0[j] += alpha * ph.d0[j] + omega * sh.d0[j];
      r[j] = s[j] - omega * t[j];
    }
    /* are we done? */
    *normres = norm2SF (r, ntotal);
    if (output == 1) {
      printf ("bicgstab: %5d  %10.3e  %10.3e  %10.3e  %10.3e\n",
	      ii + 1, (double) *normres, (double)alpha, (double)beta, (double)omega);
      fflush(stdout);
    }
    if (*normres <= tol)
      break;
    rho1 = rho;
    if (fabs (omega) < omegatol)
      break;

  }

  /* free temporary storage */
  free_dvectorSF (r, 0, ntotal - 1);
  free_dvectorSF (p, 0, ntotal - 1);
/*      free_dvectorSF(ph,  0, ntotal-1);*/
  free_derivs (&ph, ntotal);
  free_dvectorSF (rt, 0, ntotal - 1);
  free_dvectorSF (s, 0, ntotal - 1);
/*      free_dvectorSF(sh,  0, ntotal-1);*/
  free_derivs (&sh, ntotal);
  free_dvectorSF (t, 0, ntotal - 1);
  free_dvectorSF (vv, 0, ntotal - 1);

  free_dvectorSF (F, 0, ntotal - 1);
  free_derivs (&u, ntotal);

  free_dmatrixSF (JFD, 0, ntotal - 1, 0, maxcol - 1);
  free_imatrixSF (cols, 0, ntotal - 1, 0, maxcol - 1);
  free_ivectorSF (ncols, 0, ntotal - 1);

  /* iteration failed */
  if (ii > itmax)
    return -1;

  /* breakdown */
  if (fabs (rho) < rhotol)
    return -10;
  if (fabs (omega) < omegatol)
    return -11;

  /* success! */
  return ii + 1;
}

/* -------------------------------------------------------------------*/
void
NewtonSF (CCTK_POINTER_TO_CONST const cctkGH,
        int const nvar, int const n1, int const n2, int const n3,
	derivs v,
        CCTK_REAL const tol, int const itmax)
{
  DECLARE_CCTK_PARAMETERS;

  int ntotal = n1 * n2 * n3 * nvar, ii, it;
  CCTK_REAL *F, dmax, normres;
  derivs u, dv;

  F = dvectorSF (0, ntotal - 1);
  allocate_derivs (&dv, ntotal);
  allocate_derivs (&u, ntotal);

/*         TestRelax(nvar, n1, n2, n3, v, dv.d0); */
  it = 0;
  dmax = 1;
  while (dmax > tol && it < itmax)
  {
    CCTK_VInfo (CCTK_THORNSTRING, "Attempting to find bare masses. %d", it);
    if (it == 0)
    {
      F_of_vSF (cctkGH, nvar, n1, n2, n3, v, F, u);
      dmax = norm_infSF (F, ntotal);
    }
#pragma omp parallel for
    for (int j = 0; j < ntotal; j++)
      dv.d0[j] = 0;

    if(verbose==1){
      printf ("Newton: it=%d \t |F|=%e\n", it, (double)dmax);
      printf ("bare mass: mp=%g \t mm=%g\n", (double) par_m_plus, (double) par_m_minus);
    }

    fflush(stdout);
    ii =
      bicgstabSF (cctkGH,
                nvar, n1, n2, n3, v, dv, verbose, 100, dmax * 1.e-3, &normres);
#pragma omp parallel for
    for (int j = 0; j < ntotal; j++)
      v.d0[j] -= dv.d0[j];
    F_of_vSF (cctkGH, nvar, n1, n2, n3, v, F, u);
    dmax = norm_infSF (F, ntotal);
    it += 1;
  }
  if (itmax==0)
  {
      F_of_vSF (cctkGH, nvar, n1, n2, n3, v, F, u);
      dmax = norm_infSF (F, ntotal);
  }

  //if(verbose==1)
  if(1)
    printf ("Newton: it=%d \t |F|=%e \n", it, (double)dmax);

  fflush(stdout);

  free_dvectorSF (F, 0, ntotal - 1);
  free_derivs (&dv, ntotal);
  free_derivs (&u, ntotal);
}

/* -------------------------------------------------------------------*/
