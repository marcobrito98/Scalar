#define StencilSize 19
#define N_PlaneRelax 1
#define NRELAX 200
#define Step_Relax 1

typedef struct DERIVS {
  CCTK_REAL *d0, *d1, *d2, *d3, *d11, *d12, *d13, *d22, *d23, *d33;
} derivs;

/* Routines in  "TwoPunctures_BBHSF.c"*/
CCTK_REAL TestSolution(CCTK_REAL A, CCTK_REAL B, CCTK_REAL X, CCTK_REAL R,
                       CCTK_REAL phi);
void TestVector_w(CCTK_REAL *par, int nvar, int n1, int n2, int n3,
                  CCTK_REAL *w);

/* Routines in  "FuncAndJacobian.c"*/
int BBHSF_Index(int ivar, int i, int j, int k, int nvar, int n1, int n2,
                int n3);
void BBHSF_allocate_derivs(derivs *v, int n);
void BBHSF_free_derivs(derivs *v, int n);
void BBHSF_Derivatives_AB3(int nvar, int n1, int n2, int n3, derivs v);
void BBHSF_F_of_v(CCTK_POINTER_TO_CONST cctkGH, int nvar, int n1, int n2,
                  int n3, derivs v, CCTK_REAL *F, derivs u);
void BBHSF_J_times_dv(int nvar, int n1, int n2, int n3, derivs dv,
                      CCTK_REAL *Jdv, derivs u);
void BBHSF_JFD_times_dv(int i, int j, int k, int nvar, int n1, int n2, int n3,
                        derivs dv, derivs u, CCTK_REAL *values);
void BBHSF_SetMatrix_JFD(int nvar, int n1, int n2, int n3, derivs u, int *ncols,
                         int **cols, CCTK_REAL **Matrix);
CCTK_REAL BBHSF_PunctEvalAtArbitPosition(CCTK_REAL *v, int ivar, CCTK_REAL A,
                                         CCTK_REAL B, CCTK_REAL phi, int nvar,
                                         int n1, int n2, int n3);
void BBHSF_calculate_derivs(int i, int j, int k, int ivar, int nvar, int n1,
                            int n2, int n3, derivs v, derivs vv);
CCTK_REAL BBHSF_interpol(CCTK_REAL a, CCTK_REAL b, CCTK_REAL c, derivs v);
CCTK_REAL BBHSF_PunctTaylorExpandAtArbitPosition(int ivar, int nvar, int n1,
                                                 int n2, int n3, derivs v,
                                                 CCTK_REAL x, CCTK_REAL y,
                                                 CCTK_REAL z);
CCTK_REAL BBHSF_PunctIntPolAtArbitPosition(int ivar, int nvar, int n1, int n2,
                                           int n3, derivs v, CCTK_REAL x,
                                           CCTK_REAL y, CCTK_REAL z);
void BBHSF_SpecCoef(int n1, int n2, int n3, int ivar, CCTK_REAL *v,
                    CCTK_REAL *cf);
CCTK_REAL PunctEvalAtArbitPositionFast(CCTK_REAL *v, int ivar, CCTK_REAL A,
                                       CCTK_REAL B, CCTK_REAL phi, int nvar,
                                       int n1, int n2, int n3);
CCTK_REAL BBHSF_PunctIntPolAtArbitPositionFast(int ivar, int nvar, int n1,
                                               int n2, int n3, derivs v,
                                               CCTK_REAL x, CCTK_REAL y,
                                               CCTK_REAL z);

/* Routines in  "CoordTransf.c"*/
void BBHSF_AB_To_XR(int nvar, CCTK_REAL A, CCTK_REAL B, CCTK_REAL *X,
                    CCTK_REAL *R, derivs U);
void BBHSF_C_To_c(int nvar, CCTK_REAL X, CCTK_REAL R, CCTK_REAL *x,
                  CCTK_REAL *r, derivs U);
void BBHSF_rx3_To_xyz(int nvar, CCTK_REAL x, CCTK_REAL r, CCTK_REAL phi,
                      CCTK_REAL *y, CCTK_REAL *z, derivs U);

/* Routines in  "Equations.c"*/
CCTK_REAL BBHSF_BY_KKofxyz(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z);
void BBHSF_BY_Aijofxyz(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
                       CCTK_REAL Aij[3][3]);
void BBHSF_NonLinEquations(CCTK_REAL rho_adm, CCTK_REAL A, CCTK_REAL B,
                           CCTK_REAL X, CCTK_REAL R, CCTK_REAL x, CCTK_REAL r,
                           CCTK_REAL phi, CCTK_REAL y, CCTK_REAL z, derivs U,
                           CCTK_REAL *values);
void BBHSF_LinEquations(CCTK_REAL A, CCTK_REAL B, CCTK_REAL X, CCTK_REAL R,
                        CCTK_REAL x, CCTK_REAL r, CCTK_REAL phi, CCTK_REAL y,
                        CCTK_REAL z, derivs dU, derivs U, CCTK_REAL *values);

/* Routines in  "Newton.c"*/
void BBHSF_TestRelax(CCTK_POINTER_TO_CONST cctkGH, int nvar, int n1, int n2,
                     int n3, derivs v, CCTK_REAL *dv);
void BBHSF_Newton(CCTK_POINTER_TO_CONST cctkGH, int nvar, int n1, int n2,
                  int n3, derivs v, CCTK_REAL tol, int itmax);

/* Routines in "SF_density.c"*/
void SF_Initialize(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL *phi_re,
                   CCTK_REAL *phi_im);
void conf_flat_analytic_SF_source_term(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
                                       CCTK_REAL phi_re, CCTK_REAL phi_im,
                                       CCTK_REAL *sourceSF1,
                                       CCTK_REAL *sourceSF2,
                                       CCTK_REAL Dphi_re_dx[3],
                                       CCTK_REAL Dphi_im_dx[3]);
