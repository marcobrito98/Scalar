/* TwoPuncturesSF:  File  "utilities.h"*/

#include <math.h>

#include "cctk.h"

#define Pi  3.14159265358979323846264338328
#define Pih 1.57079632679489661923132169164	/* Pi/2*/
#define Piq 0.78539816339744830961566084582	/* Pi/4*/

#define TINY 1.0e-20
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

#define nrerror TP_nrerror
#define ivectorSF TP_ivectorSF
#define dvectorSF TP_dvectorSF
#define imatrixSF TP_imatrixSF
#define dmatrixSF TP_dmatrixSF
#define d3tensorSF TP_d3tensorSF
#define free_ivectorSF TP_free_ivectorSF
#define free_dvectorSF TP_free_dvectorSF
#define free_imatrixSF TP_free_imatrixSF
#define free_dmatrixSF TP_free_dmatrixSF
#define free_d3tensorSF TP_free_d3tensorSF

#define minimum2SF TP_minimum2SF
#define minimum3SF TP_minimum3SF
#define maximum2SF TP_maximum2SF
#define maximum3SF TP_maximum3SF
#define pow_intSF TP_pow_intSF

#define chebft_ZerosSF TP_chebft_ZerosSF
#define chebft_ExtremesSF TP_chebft_ExtremesSF
#define chderSF TP_chderSF
#define chebevSF TP_chebevSF
#define fourftSF TP_fourftSF
#define fourderSF TP_fourderSF
#define fourder2SF TP_fourder2SF
#define fourevSF TP_fourevSF

#define norm1SF TP_norm1SF
#define norm2SF TP_norm2SF
#define scalarproductSF TP_scalarproductSF

void nrerror (char error_text[]);
int *ivectorSF (long nl, long nh);
CCTK_REAL *dvectorSF (long nl, long nh);
int **imatrixSF (long nrl, long nrh, long ncl, long nch);
CCTK_REAL **dmatrixSF (long nrl, long nrh, long ncl, long nch);
CCTK_REAL ***d3tensorSF (long nrl, long nrh, long ncl, long nch, long ndl,
		    long ndh);
void free_ivectorSF (int *v, long nl, long nh);
void free_dvectorSF (CCTK_REAL *v, long nl, long nh);
void free_imatrixSF (int **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrixSF (CCTK_REAL **m, long nrl, long nrh, long ncl, long nch);
void free_d3tensorSF (CCTK_REAL ***t, long nrl, long nrh, long ncl, long nch,
		    long ndl, long ndh);

int minimum2SF (int i, int j);
int minimum3SF (int i, int j, int k);
int maximum2SF (int i, int j);
int maximum3SF (int i, int j, int k);
int pow_intSF (int mantisse, int exponent);

void chebft_ZerosSF (CCTK_REAL u[], int n, int inv);
void chebft_ExtremesSF (CCTK_REAL u[], int n, int inv);
void chderSF (CCTK_REAL *c, CCTK_REAL *cder, int n);
CCTK_REAL chebevSF (CCTK_REAL a, CCTK_REAL b, CCTK_REAL c[], int m, CCTK_REAL x);
void fourftSF (CCTK_REAL *u, int N, int inv);
void fourderSF (CCTK_REAL u[], CCTK_REAL du[], int N);
void fourder2SF (CCTK_REAL u[], CCTK_REAL d2u[], int N);
CCTK_REAL fourevSF (CCTK_REAL *u, int N, CCTK_REAL x);


CCTK_REAL norm1SF (CCTK_REAL *v, int n);
CCTK_REAL norm2SF (CCTK_REAL *v, int n);
CCTK_REAL scalarproductSF (CCTK_REAL *v, CCTK_REAL *w, int n);
