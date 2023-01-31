/* TwoPuncturesSF:  File  "TwoPuncturesSF.h"*/

#define StencilSize 19
#define N_PlaneRelax 1
#define NRELAX 200
#define Step_Relax 1

typedef struct DERIVS
{
  CCTK_REAL *d0, *d1, *d2, *d3, *d11, *d12, *d13, *d22, *d23, *d33;
} derivs;

/*
Files of "TwoPuncturesSF":
	TwoPuncturesSF.c
	FuncAndJacobian.c
	CoordTransf.c
	Equations.c
	Newton.c
	utilities.c (see utilities.h)
**************************
*/

/* Routines in  "TwoPuncturesSF.c"*/
CCTK_REAL TestSolutionSF (CCTK_REAL A, CCTK_REAL B, CCTK_REAL X, CCTK_REAL R, CCTK_REAL phi);
void TestVector_wSF (CCTK_REAL *par, int nvar, int n1, int n2, int n3, CCTK_REAL *w);

/* Routines in  "FuncAndJacobian.c"*/
int IndexSF (int ivar, int i, int j, int k, int nvar, int n1, int n2, int n3);
void allocate_derivsSF (derivs * v, int n);
void free_derivsSF (derivs * v, int n);
void Derivatives_AB3SF (int nvar, int n1, int n2, int n3, derivs v);
void F_of_vSF (CCTK_POINTER_TO_CONST cctkGH,
       int nvar, int n1, int n2, int n3, derivs v,
	     CCTK_REAL *F, derivs u);
void J_times_dvSF (int nvar, int n1, int n2, int n3, derivs dv,
		 CCTK_REAL *Jdv, derivs u);
void JFD_times_dvSF (int i, int j, int k, int nvar, int n1,
		   int n2, int n3, derivs dv, derivs u, CCTK_REAL *values);
void SetMatrix_JFDSF (int nvar, int n1, int n2, int n3,
		    derivs u, int *ncols, int **cols, CCTK_REAL **Matrix);
CCTK_REAL PunctEvalAtArbitPositionSF (CCTK_REAL *v, int ivar, CCTK_REAL A, CCTK_REAL B, CCTK_REAL phi,
				 int nvar, int n1, int n2, int n3);
void calculate_derivsSF (int i, int j, int k, int ivar, int nvar, int n1,
		       int n2, int n3, derivs v, derivs vv);
CCTK_REAL interpolSF (CCTK_REAL a, CCTK_REAL b, CCTK_REAL c, derivs v);
CCTK_REAL PunctTaylorExpandAtArbitPositionSF (int ivar, int nvar, int n1,
                                         int n2, int n3, derivs v, CCTK_REAL x,
                                         CCTK_REAL y, CCTK_REAL z);
CCTK_REAL PunctIntPolAtArbitPositionSF (int ivar, int nvar, int n1,
				   int n2, int n3, derivs v, CCTK_REAL x,
				   CCTK_REAL y, CCTK_REAL z);
void SpecCoefSF(int n1, int n2, int n3, int ivar, CCTK_REAL *v, CCTK_REAL *cf);
CCTK_REAL PunctEvalAtArbitPositionFastSF (CCTK_REAL *v, int ivar, CCTK_REAL A, CCTK_REAL B, CCTK_REAL phi,
                                 int nvar, int n1, int n2, int n3);
CCTK_REAL PunctIntPolAtArbitPositionFastSF (int ivar, int nvar, int n1,
                                   int n2, int n3, derivs v, CCTK_REAL x,
                                   CCTK_REAL y, CCTK_REAL z);


/* Routines in  "CoordTransf.c"*/
void AB_To_XRSF (int nvar, CCTK_REAL A, CCTK_REAL B, CCTK_REAL *X,
	       CCTK_REAL *R, derivs U);
void C_To_cSF (int nvar, CCTK_REAL X, CCTK_REAL R, CCTK_REAL *x,
	     CCTK_REAL *r, derivs U);
void rx3_To_xyzSF (int nvar, CCTK_REAL x, CCTK_REAL r, CCTK_REAL phi, CCTK_REAL *y,
		 CCTK_REAL *z, derivs U);

/* Routines in  "Equations.c"*/
CCTK_REAL BY_KKofxyzSF (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z);
void BY_AijofxyzSF (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL Aij[3][3]);
void NonLinEquationsSF (CCTK_REAL rho_adm,
          CCTK_REAL A, CCTK_REAL B, CCTK_REAL X, CCTK_REAL R,
		      CCTK_REAL x, CCTK_REAL r, CCTK_REAL phi,
		      CCTK_REAL y, CCTK_REAL z, derivs U, CCTK_REAL *values);
void LinEquationsSF (CCTK_REAL A, CCTK_REAL B, CCTK_REAL X, CCTK_REAL R,
		   CCTK_REAL x, CCTK_REAL r, CCTK_REAL phi,
		   CCTK_REAL y, CCTK_REAL z, derivs dU, derivs U, CCTK_REAL *values);

/* Routines in  "Newton.c"*/
void TestRelaxSF (CCTK_POINTER_TO_CONST cctkGH,
                int nvar, int n1, int n2, int n3, derivs v, CCTK_REAL *dv);
void NewtonSF (CCTK_POINTER_TO_CONST cctkGH,
             int nvar, int n1, int n2, int n3, derivs v,
	           CCTK_REAL tol, int itmax);

/* Routines in  "ScalarSource.c"*/
void SF_Gaussian(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL *phi2, CCTK_REAL *dphi2);
CCTK_REAL LinSrcSF (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL psi);
CCTK_REAL NonLinSrcSF (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL psi);

/* 
 27: -1.325691774825335e-03
 37: -1.325691778944117e-03
 47: -1.325691778942711e-03
 
 17: -1.510625972641537e-03
 21: -1.511443006977708e-03
 27: -1.511440785153687e-03
 37: -1.511440809549005e-03
 39: -1.511440809597588e-03
 */
