/* intial data thorn: ID_SF_BS */
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
void ID_SF_BS(CCTK_ARGUMENTS);
void
ID_SF_BS (CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

//  CCTK_INFO("=== Begin ID_SF_BS initial data ===");

  /*=== define parameters ====================================*/
  // Scalar field mass parameter
  CCTK_REAL mu2, mu4;
  mu2 = mu  * mu;
  mu4 = mu2 * mu2;

  //Note: we use dimensionless spin parameter
  CCTK_REAL spin, spin2;
  spin  = spin_plus / m_plus;
  spin2 = spin*spin;

  // horizon radius
  CCTK_REAL BB, rBLp, rBLm;
  BB   = sqrt( 1.0 - spin2 );
  rBLp = m_plus + BB;
  rBLm = m_plus - BB;

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
     // scalar field vars
     CCTK_REAL psi_re, psi_im, psit_re, psit_im;
     psi_re  = 0.0;
     psi_im  = 0.0;
     psit_re = 0.0;
     psit_im = 0.0;

     // radial and theta functions
     CCTK_REAL RR11, RI11, SR11, SI11;
     RR11 = 0.0;
     RI11 = 0.0;
     SR11 = 0.0;
     SI11 = 0.0;
     /*========================================================*/

     /*========================================================*/

     /*=== define radius, r_{BL} and angles ===================*/
     // positions
     CCTK_REAL xp[3];
//     xp[0] = x[ind];
//     xp[1] = y[ind];
//     xp[2] = z[ind];

     xp[0] = x[ind] - pos_plus[0]; 
     xp[1] = y[ind] - pos_plus[1];
     xp[2] = z[ind] - pos_plus[2];

     // coordinate radius and polar radial coordinate
     CCTK_REAL rr, rr2;
     rr2 = xp[0] * xp[0] + xp[1] * xp[1] + xp[2] * xp[2];
     if( rr2 < eps_r2 ) rr2 = eps_r2;
     rr = sqrt( rr2 );

     CCTK_REAL rho, rho2;
     rho2 =  xp[0]*xp[0] + xp[1]*xp[1];
     if( rho2 < eps_r2 ) rho2 = eps_r2;
     rho  = sqrt( rho2 );

     // Boyer-Lindquist radial coordinate
     CCTK_REAL rBL;
     rBL  = rr * ( 1.0 + 0.25 * rBLp / rr ) * ( 1.0 + 0.25 * rBLp / rr );
     if( rBL < eps_r ) rBL = eps_r;

     // differences rBL - rBL_{\pm}
     CCTK_REAL Rp, Rm;

     Rp = rBL - rBLp;
     if( abs( Rp ) < eps_r ) Rp = eps_r;

     Rm = rBL - rBLm;
     if( abs(Rm) < eps_r ) Rm = eps_r;

     // angles
     CCTK_REAL CT, CT2, CP;
     CCTK_REAL SP;

     CT  = xp[2] / rr;
     CT2 = CT*CT;
     CP  = xp[0] / rho;

     SP  = xp[1] / rho;
     /*========================================================*/

     /*=== define coefficients, eqs (34) in Dolan '07 =========*/
     // frequency vars
     CCTK_REAL wcrit, wR2, wR3, wR4, wI2, wI3, wI4;
     wcrit = spin / ( 2*rBLp );
     wR2   = wR  * wR;
     wR3   = wR2 * wR;
     wR4   = wR2 * wR2;
     wI2   = wI  * wI;
     wI3   = wI2 * wI;
     wI4   = wI2 * wI2;


     CCTK_REAL sigmaR, sigmaI;
     sigmaR = ( BB + 1.0 ) * ( wR - wcrit ) / BB;
     sigmaI = ( BB + 1.0 ) * wI		    / BB; 


     CCTK_REAL qR, qI;
     // NOTE: this uses a series expansion around $\omega_{I} \sim 0$ 

     if( wR2 < mu2 )
     {	
       qR = - ( mu4 + wR4 + 0.5 * mu2 * ( wI2 - 4.0 * wR2 ) ) / sqrt( pow( mu2 - wR2, 3 ) );
       qI =   wR * wI / sqrt( mu2 - wR2 );

     }
     else if( mu2 < wR2 )
     {
       qR = wR * wI / sqrt( wR2 - mu2 );
       qI = - ( mu4 + wR4 + 0.5 * mu2 * ( wI2 - 4.0 * wR2 ) ) / sqrt( pow( wR2 - mu2, 3 ) );

     }
     else if( mu2 == wR2 && wI > 0 )
     {
       qR = - sqrt( wI / mu ) * ( mu + 0.25*wI );
       qI =   sqrt( wI / mu ) * ( mu - 0.25*wI );

     }
     else if( mu2 == wR2 && wI < 0 )
     {
       qR = sqrt( abs(wI) / mu ) * ( mu - 0.25*wI );
       qI = sqrt( abs(wI) / mu ) * ( mu + 0.25*wI );

     }

     CCTK_REAL qR2, qI2, q2, q4;
     qR2 = qR * qR;
     qI2 = qI * qI;
     q2  = qR2 + qI2;
     q4  = q2 * q2;


     CCTK_REAL chiR, chiI;
     chiR = - 4.0*qI * wR * wI / q2 + qR * ( mu2 + 2.0*( wI2 - wR2 ) ) / q2;
     chiI = - 4.0*qR * wR * wI / q2 - qI * ( mu2 + 2.0*( wI2 - wR2 ) ) / q2;

     /*========================================================*/


     /*=== X1R, X1I with X1 = Rp^{-i sigma } ==================*/
     CCTK_REAL X1R, X1I;

     if( Rp > 0.0 )
     {
       X1R =   pow( Rp, sigmaI ) * cos( sigmaR * log( Rp ) );
       X1I = - pow( Rp, sigmaI ) * sin( sigmaR * log( Rp ) ); 
     }
     else if( Rp < 0.0 )
     {
       X1R = pow( abs(Rp), sigmaI ) * exp( Pi*sigmaR ) * cos( Pi*sigmaI - sigmaR * log( abs(Rp) ) );
       X1I = pow( abs(Rp), sigmaI ) * exp( Pi*sigmaR ) * sin( Pi*sigmaI - sigmaR * log( abs(Rp) ) );
     }
     else
     {
       X1R = 0.0;
       X1I = 0.0;
     }

     /*========================================================*/

     /*=== X2R, X2I with X2 = Rm^{i sigma + chi -1} ===========*/
     CCTK_REAL X2R, X2I;
     CCTK_REAL fac;
     fac = chiR - sigmaI - 1.0;

     if( Rm > 0.0 )
     {
       X2R =   pow( Rm, fac ) * cos( (sigmaR + chiI) * log(Rm) ); 
       X2I =   pow( Rm, fac ) * sin( (sigmaR + chiI) * log(Rm) );
     }
     else if( Rm < 0.0 )
     {
       X2R = - pow( abs(Rm), fac ) * exp( - Pi * ( chiI + sigmaR ) ) 
		* cos( Pi * ( chiR - sigmaI ) + ( chiI + sigmaR ) * log(abs(Rm)) );
       X2I =   pow( abs(Rm), fac ) * exp( - Pi * ( chiI + sigmaR ) )
		* sin( Pi * ( sigmaI - chiR ) - ( chiI + sigmaR ) * log(abs(Rm)) );
     }
     else
     {
       X2R = 0.0;
       X2I = 0.0;
     }

     /*========================================================*/

     /*=== X3R, X3I with X3 = exp(q rBL) ======================*/
     CCTK_REAL X3R, X3I;

     X3R = exp( qR*rBL ) * cos( qI*rBL );
     X3I = exp( qR*rBL ) * sin( qI*rBL );

     /*========================================================*/

     /*=== X4R, X4I with X4 = Sum_n An (Rp/Rm)^n ==============*/
     CCTK_REAL X4R, X4I;

	/*--- Lam11, eigenvalue of s=0 spheroidal harmonic ----*/

	CCTK_REAL ff[4];

	ff[0] =   2.0;          // f_{0}
	ff[1] = - 1.0 / 5;      // f_{2}
	ff[2] = - 4.0 / 875;    // f_{4}
	ff[3] = - 8.0 / 65625;  // f_{6}

	CCTK_REAL CCR, CCR2, CCR3, CCI, CCI2;
        CCR = - spin2 * ( mu2 + wI2 - wR2 );
        CCI = 2.0 * spin2 *  wI * wR;

	CCR2 = CCR  * CCR;
	CCR3 = CCR2 * CCR;
 	CCI2 = CCI  * CCI;


	CCTK_REAL LamR11, LamI11;
	LamR11 = ff[0] + ff[2] * ( CCR2 - CCI2 ) + CCR3 * ff[3] + CCR * ( ff[1] - 3 * CCI2 * ff[3] );
	LamI11 = CCI * ( ff[1] + 2*CCR * ff[2] + ff[3] * ( 3*CCR2 - CCI2 ) );

        /*-----------------------------------------------------*/

        /*--- C0, C1, C2, C3, C4, Eqs. 40 - 44 in Dolan '07 ---*/

	CCTK_REAL C0R, C0I;

	C0R = 1.0 + 2.0 * wI * ( 1.0 + BB ) / BB;
	C0I = ( spin - 2.0 * wR * ( 1.0 + BB ) ) / BB;

	/*-----------------------------------------------------*/
	CCTK_REAL C1R, C1I;

	C1R = - 4.0 + 2*qR * ( 1.0 + 2.0*BB ) - 4*wI * ( 1.0 + 1.0/BB )
              + 2*qR * ( wI2 - wR2 ) / q2 - 4*qI*wI*wR / q2;

	C1I =   2*qI * ( 1.0 + 2.0*BB ) + 4*wR * ( 1.0 + 1.0/BB ) - 2*spin / BB
              + 2*qI * ( wR2 - wI2 ) / q2 - 4*qR*wI*wR / q2;

	/*-----------------------------------------------------*/
	CCTK_REAL C2R, C2I;

	C2R =   3.0 - 2*qR + 2*wI * ( 1.0 + 1.0/BB )
              + 2.0*qR * ( wR2 - wI2 ) / q2 + 4*qI*wI*wR / q2;

	C2I = - 2.0*qI - 2*wR * ( 1.0 + 1.0/BB ) + spin/BB
              + 2.0*qI * ( wI2 - wR2 ) / q2 + 4*qR*wI*wR / q2;

	/*-----------------------------------------------------*/
	CCTK_REAL C3R, C3I;

	C3R = - LamR11 - 1.0
              + qR + 6 * ( wR2 - wI2 + qI*wR + qR*wI ) - 2*qI*spin - 2*wI
              + q2 * ( 2.0 - spin2 ) + 2*qR2 * ( spin2 - 2.0 )
              + 2*BB * ( wR2 - wI2 ) + 2*BB * ( q2 + qR ) + 4*BB * ( qR*wI + qI*wR - qR2 )
              - qI*spin / BB + 4 * ( wR2 - wI2 ) / BB + 2 * ( qR*wI + qI*wR - spin*wR - wI ) / BB
              + qR * ( wI2 - wR2 ) / q2  - 6*wI*wR * ( wI*qI + wR*qR ) / q2
              + 2 * ( qR*wI3 + qI*wR3 - qI*wI*wR ) / q2
              + qI*spin * ( wI2 - wR2 ) / (BB*q2) - 6*wI*wR * ( wI*qI + wR*qR ) / (BB*q2)
              + 2 * ( qR*wI3 + qI*wR3 + qR*spin*wI*wR ) / (BB*q2);

	C3I = - LamI11
              +  qI * ( 1.0 - 4*qR ) + 2*qR*spin * ( 1.0 + qI*spin )
              + 6 * ( qI*wI - qR*wR ) + 2*wR * ( 1.0 + 6*wI )
              + 2*BB*qI * ( 1.0 - 2*qR ) + 4*BB * ( qI*wI - qR*wR + wI*wR )
              + spin * ( qR - 1.0 - 2*wI ) / BB
              + 2 * ( qI*wI - qR*wR ) / BB + 2*wR * ( 1.0 + 4*wI ) / BB
              + qI * ( wR2 - wI2 ) / q2 + 6*wI*wR * ( wR*qI - wI*qR ) / q2
              + 2 * ( qR*wR3 - qI*wI3 - qR*wI*wR ) / q2
              + qR*spin * ( wI2 - wR2 ) / (BB*q2) + 6*wI*wR * ( wR*qI - wI*qR ) / (BB*q2)
              + 2 * ( qR*wR3 - qI*wI3 - qI*spin*wI*wR ) / (BB*q2);

	/*-----------------------------------------------------*/
	CCTK_REAL C4R, C4I;

	C4R = - q2 + 2 * ( qR2 - qR*wI + wI2 - wR2 - qI*wR )
              + spin * ( qI + 2*wR ) / BB
              + 4 * ( wI2 - wR2 ) / BB - 2 * ( qI*wR + qR*wI )/ BB
              - ( wR4 + wI4 ) / q2 - 2 * ( qR*wI3 + qI*wR3 ) / q2
              + 6*wI*wR * ( wI*qI + wR*qR + wI*wR ) / q2
              + qI*spin * ( wR2 - wI2 ) / (BB*q2)
              + 6*wI*wR * ( wI*qI + wR*qR ) / (BB*q2)
              - 2 * ( qR*wI3 + qI*wR3 + qR*spin*wI*wR ) / (BB*q2)
              + 2*qR2 * ( wI4 + wR4 - 6*wI2*wR2 ) / q4
              + 8*qI*qR*wI*wR * ( wR2 - wI2 ) / q4;

	C4I =   2 * ( qI*qR + qR*wR - qI*wI - 2*wI*wR )
              + spin * ( 2*wI - qR ) / BB + 2 * ( qR*wR - qI*wI - 4*wI*wR ) / BB
              + 2*wI3 * ( qI + 2*wR ) / q2 - 2*wR3 * ( qR + 2*wI ) / q2
              + 6*wI*wR * ( wI*qR - wR*qI ) / q2
              + qR*spin * ( wR2 - wI2 ) / (BB*q2)
              + 6*wI*wR * ( wI*qR - wR*qI ) / (BB*q2)
              + 2 * ( qI*wI3 + qI*spin*wI*wR - qR*wR3 ) / (BB*q2)
              - 2*qI*qR * ( wI4 + wR4 ) / q4 + 12*qI*qR*wI2*wR2 / q4
              + 8*qR2*wI*wR * ( wR2 - wI2 ) / q4;

	/*-----------------------------------------------------*/

	/*--- alpha, beta, gamma, Eqs. 37 - 39 in Dolan '07 ---*/
	// alpha
	CCTK_REAL alphR[n_sum], alphI[n_sum], alph2[n_sum];
	for( int m = 0; m < n_sum; ++m )
	{

	  alphR[m] = ( C0R + m ) * ( m + 1 );
          alphI[m] = C0I	 * ( m + 1 );
	  alph2[m] = alphR[m]*alphR[m] + alphI[m]*alphI[m];

	}

	/*-----------------------------------------------------*/
	// beta
	CCTK_REAL betaR[n_sum], betaI[n_sum];
	for( int m = 0; m < n_sum; ++m )
	{

          betaR[m] = C3R + m * ( C1R + 2*(1-m) );
          betaI[m] = C3I + m * C1I;

	}

	/*-----------------------------------------------------*/
	// gamma
    	CCTK_REAL gamR[n_sum], gamI[n_sum];
	for( int m = 0; m < n_sum; ++m )
	{

          gamR[m] = C4R + m * ( C2R + m - 3 );
          gamI[m] = C4I + m * C2I;

	}

	/*-----------------------------------------------------*/

	/*--- AnR, AnI, Eqs. 35,36 in Dolan '07 ---------------*/
	CCTK_REAL AnR[n_sum], AnI[n_sum];

	AnR[0] = A0;
	AnI[0] = 0.0;

	AnR[1] = - AnR[0] * ( alphR[0] * betaR[0] + alphI[0] * betaI[0] ) / alph2[0];
	AnI[1] = - AnR[0] * ( alphR[0] * betaI[0] - alphI[0] * betaR[0] ) / alph2[0];

	for( int m = 2; m < n_sum; ++m)
	{

	   AnR[m] = ( - AnR[m-1] * ( alphR[m-1] * betaR[m-1] + alphI[m-1] * betaI[m-1] )
		      + AnI[m-1] * ( alphR[m-1] * betaI[m-1] - alphI[m-1] * betaR[m-1] )
		      - AnR[m-2] * ( alphR[m-1] * gamR[m-1]  + alphI[m-1] * gamI[m-1]  )
		      + AnI[m-2] * ( alphR[m-1] * gamI[m-1]  - alphI[m-1] * gamR[m-1]  )
		    ) / alph2[m-1];

	   AnI[m] = ( - AnR[m-1] * ( alphR[m-1] * betaI[m-1] - alphI[m-1] * betaR[m-1] )
		      - AnI[m-1] * ( alphR[m-1] * betaR[m-1] + alphI[m-1] * betaI[m-1] )
		      - AnR[m-2] * ( alphR[m-1] * gamI[m-1]  - alphI[m-1] * gamR[m-1]  )
		      - AnI[m-2] * ( alphR[m-1] * gamR[m-1]  + alphI[m-1] * gamI[m-1]  ) 
		    ) / alph2[m-1];

	}

	/*-----------------------------------------------------*/
     // X4R, X4I 
     X4R = 0.0;
     X4I = 0.0;

     for( int m = 0; m < n_sum; ++m )
     {

       X4R = X4R + AnR[m] * pow( Rp / Rm, m );
       X4I = X4I + AnI[m] * pow( Rp / Rm, m );

     }
     /*========================================================*/

     /*=== initialize R_{11} spherical harmonics ( in cartesian coordinates ) ===*/

     RR11 =  X1I * ( X2I*X3I*X4I - X2R*X3R*X4I - X2R*X3I*X4R - X2I*X3R*X4R )
           - X1R * ( X2R*X3I*X4I + X2I*X3R*X4I + X2I*X3I*X4R - X2R*X3R*X4R );

     RI11 =  X1R * (-X2I*X3I*X4I + X2R*X3R*X4I + X2R*X3I*X4R + X2I*X3R*X4R )
           - X1I * ( X2R*X3I*X4I + X2I*X3R*X4I + X2I*X3I*X4R - X2R*X3R*X4R );

     /*========================================================*/

     /*=== s=0 spheroidal harmonics ===========================*/
     // Note: see Eqs. (2.2) - (2.9) in gr-qc/0511111
	/*-----------------------------------------------------*/
	CCTK_REAL cR, cI;
	if( wR2 < mu2 )
	{
	  cR = spin * wI * wR / sqrt( mu2 - wR2 );
	  cI = spin * ( mu4 + wR4 + 0.5 * mu2 * ( wI2 - 4.0*wR2 ) ) / sqrt( pow( mu2 - wR2, 3 ) );
	}
	else if( mu2 < wR2 )
	{
	  cR = spin * ( mu4 + wR4 + 0.5 * mu2 * ( wI2 - 4.0*wR2 ) ) / sqrt( pow( wR2 - mu2, 3 ) );
	  cI = spin * wI * wR / sqrt( wR2 - mu2 );
	}
     	else if( mu2 == wR2 && wI > 0 )
	{
	  cR =   spin * sqrt( wI / mu ) * ( mu - 0.25*wI );
	  cI =   spin * sqrt( wI / mu ) * ( mu + 0.25*wI );
	}
     	else if( mu2 == wR2 && wI < 0 )
	{
	  cR =   spin * sqrt( abs(wI) / mu ) * ( mu + 0.25*wI );
	  cI = - spin * sqrt( abs(wI) / mu ) * ( mu - 0.25*wI );
	}


	CCTK_REAL cR2, cI2;
	cR2 = cR * cR;
	cI2 = cI * cI;

	/*-----------------------------------------------------*/
	CCTK_REAL T1R, T1I;

        T1R = exp( cR * CT ) * cos( cI * CT );
        T1I = exp( cR * CT ) * sin( cI * CT );

	/*-----------------------------------------------------*/
	//Note,  we redefine alph, beta, gam here!!!
	CCTK_REAL aa[n_sum];
	CCTK_REAL bbR[n_sum], bbI[n_sum];
	CCTK_REAL ggR[n_sum], ggI[n_sum];

	for( int p = 0; p < n_sum; ++p )
	{

          aa[p]  = -2.0 * ( p+1 ) * ( p+2 );
          bbR[p] = 2.0 + cI2 - 4 * cR - cR2 + 3 * p - 4 * cR * p + p*p - LamR11;
          bbI[p] = -4 * cI - 2 * cI * cR - 4 * cI * p - LamI11;
          ggR[p] = 2.0 * (p+1) * cR;
          ggI[p] = 2.0 * (p+1) * cI;

	}

	/*-----------------------------------------------------*/
	CCTK_REAL ApR[n_sum], ApI[n_sum];

	ApR[0] = Ap0;
	ApI[0] = 0.0;

	ApR[1] = - ApR[0] * bbR[0] / aa[0];
	ApI[1] = - ApR[0] * bbI[0] / aa[0];

	for( int p = 2; p < n_sum; ++p )
	{

	   ApR[p] = ( - ApR[p-1] * bbR[p-1] + ApI[p-1] * bbI[p-1]
		      - ApR[p-2] * ggR[p-1] + ApI[p-2] * ggI[p-1]
		    ) / aa[p-1];

	   ApI[p] = ( - ApR[p-1] * bbI[p-1] - ApI[p-1] * bbR[p-1]
		      - ApR[p-2] * ggI[p-1] - ApI[p-2] * ggR[p-1]
		    ) / aa[p-1];

	}

	/*-----------------------------------------------------*/
	CCTK_REAL  T2R, T2I;
        T2R = 0.0;
        T2I = 0.0;

	for( int p = 0; p < n_sum; ++p )
	{

          T2R = T2R + ApR[p] * pow( 1 + CT, p );
          T2I = T2I + ApI[p] * pow( 1 + CT, p );

	}

	/*-----------------------------------------------------*/

     SR11 = sqrt( 1.0 - CT2 ) * ( T1R * T2R - T1I * T2I );
     SI11 = sqrt( 1.0 - CT2 ) * ( T1R * T2I + T1I * T2R );

     /*========================================================*/

     /*=== initialize Psi =====================================*/
     // Set SF to zero inside horizon to avoid divergences near origin
     CCTK_REAL rH = ( 1.0 + BB ) / 4.0;

     if( rr < rH )
     {
       psi_re = 0.0;
       psi_im = 0.0;
     } 
     else
     {
       psi_re =  ( RR11 * SR11 - RI11 * SI11 ) * CP
               - ( RR11 * SI11 + RI11 * SR11 ) * SP;

       psi_im =  ( RR11 * SI11 + RI11 * SR11 ) * CP
               + ( RR11 * SR11 - RI11 * SI11 ) * SP;
     }

     /*========================================================*/

     /*=== calc momentum Kphi =================================*/
     // implement first part of Kphi
     // Kphi = -1/(2\alpha) ( \p_{t} - \Lie_{\beta} ) \phi
     // ASSUME FOR SIMPLICITY: \alpha=1, \beta^i = 0

     psit_re = - 0.5 * ( wR * psi_im + wI * psi_re );
     psit_im =   0.5 * ( wR * psi_re - wI * psi_im );

     /*========================================================*/

     /*=== write to grid functions ============================*/

     phi1[ind]  = psi_re;
     phi2[ind]  = psi_im;

     Kphi1[ind] = psit_re;
     Kphi2[ind] = psit_im;

     /*========================================================*/

    } /* for i */
   }  /* for j */
  }   /* for k */
/*=== end of loops over grid =================================*/
/*============================================================*/


//  CCTK_INFO("=== End ID_SF_BS initial data ===");

}
/* -------------------------------------------------------------------*/
