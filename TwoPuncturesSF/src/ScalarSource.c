/* TwoPuncturesSF:  File  "ScalarSource.c"*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk_Parameters.h"
#include "TP_utilities.h"
#include "TwoPuncturesSF.h"

void conf_flat_analytic_SF_source_term(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL *Phi_re, CCTK_REAL *Phi_im, CCTK_REAL dPhi_re[3], CCTK_REAL dPhi_im[3]) {
    /* Computes analytic scalar field profiles and its derivatives */

    DECLARE_CCTK_PARAMETERS;
    CCTK_REAL rr, rr2;
    CCTK_REAL rho, rho2;

    rr2 = x*x + y*y + z*z;
    if( rr2 < TP_Tiny )
        rr2 = TP_Tiny;
    rr  = sqrt( rr2 );

    rho2 = x*x + y*y;
    if( rho2 < TP_Tiny )
        rho2 = TP_Tiny;
    rho  = sqrt( rho2 );


    /* 
     * For initial configurations assuming Gaussian radial profile 
     */

    /* Angular dependence as l=m=0 (monopole or dipole */
    CCTK_INT monopole_statement, dipole_statement;
    monopole_statement = (
            CCTK_EQUALS(initial_scalar, "ID_SF_Gaussian") &&
            CCTK_EQUALS(scalar_GaussProfile, "single_mode") &&
            l0SF == 0 && m0SF == 0 );
    dipole_statement   = (
            CCTK_EQUALS(initial_scalar, "ID_SF_Gaussian") &&
            CCTK_EQUALS(scalar_GaussProfile, "single_mode") &&
            l0SF == 1 && m0SF == 1 );
       
    if (monopole_statement) {

        /* Real and imaginary parts of the scalar field */
        *Phi_re      = sqrt(0.25/Pi) * ampSF * exp( -pow((rr - r0SF)/widthSF, 2) ); 
        *Phi_im      = 0; 

        /* Derivatives of the real and imaginary parts of the scalar field */
        dPhi_re[0]   = *Phi_re * (-2) * (rr - r0SF)/widthSF/widthSF * x / rr;
        dPhi_re[1]   = *Phi_re * (-2) * (rr - r0SF)/widthSF/widthSF * y / rr;
        dPhi_re[2]   = *Phi_re * (-2) * (rr - r0SF)/widthSF/widthSF * z / rr;
        dPhi_im[0]   = *Phi_im * (-2) * (rr - r0SF)/widthSF/widthSF * x / rr;
        dPhi_im[1]   = *Phi_im * (-2) * (rr - r0SF)/widthSF/widthSF * y / rr;
        dPhi_im[2]   = *Phi_im * (-2) * (rr - r0SF)/widthSF/widthSF * z / rr;

    } else if (dipole_statement) {

        /* Initialize radial part of the scalar field as Gaussian */
        CCTK_REAL Phi = -sqrt(0.375/Pi) * ampSF * exp( -pow((rr - r0SF)/widthSF, 2) );

        /* Derivatives of the radial function wrt (x, y, z) */
        CCTK_REAL dPhi[3];

        dPhi[0]     = Phi * (-2.) * (rr - r0SF)/widthSF/widthSF * x / rr;
        dPhi[1]     = Phi * (-2.) * (rr - r0SF)/widthSF/widthSF * y / rr;
        dPhi[2]     = Phi * (-2.) * (rr - r0SF)/widthSF/widthSF * z / rr;


        /* Angular functions and derivatives */
        CCTK_REAL sthe, dsthe[3];
        CCTK_REAL sphi, dsphi[3];
        CCTK_REAL cphi, dcphi[3];

        /* Sin and Cos (theta, phi) */
        sthe        = rho * pow(rr, -1.);
        sphi        = y * pow(rho, -1.);
        cphi        = x * pow(rho, -1.);

        /* Derivatives of Sin and Cos (theta, phi) wrt (x, y, z) */
        dsthe[0]    =   x * z * z * pow(rho, -1.) * pow(rr, -3.);
        dsthe[1]    =   y * z * z * pow(rho, -1.) * pow(rr, -3.);
        dsthe[2]    = - z * rho * pow(rr, -3.);
        dsphi[0]    = - x * y * pow(rho, -3.);
        dsphi[1]    =   x * x * pow(rho, -3.);
        dsphi[2]    = 0;
        dcphi[0]    =   y * y * pow(rho, -3.);
        dcphi[1]    = - y * x * pow(rho, -3.);
        dcphi[2]    = 0;

        /* Output */
        /* Real and imaginary parts of the scalar field */
        *Phi_re     = Phi * sthe * cphi; 
        *Phi_im     = Phi * sthe * sphi; 

        /* Derivatives of the real and imaginary parts of the scalar field */
        dPhi_re[0]  =  dPhi[0] * sthe * cphi + Phi * dsthe[0] * cphi + Phi * sthe * dcphi[0];
        dPhi_re[1]  =  dPhi[1] * sthe * cphi + Phi * dsthe[1] * cphi + Phi * sthe * dcphi[1];
        dPhi_re[2]  =  dPhi[2] * sthe * cphi + Phi * dsthe[2] * cphi + Phi * sthe * dcphi[2];

        dPhi_im[0]  =  dPhi[0] * sthe * sphi + Phi * dsthe[0] * sphi + Phi * sthe * dsphi[0];
        dPhi_im[1]  =  dPhi[1] * sthe * sphi + Phi * dsthe[1] * sphi + Phi * sthe * dsphi[1];
        dPhi_im[2]  =  dPhi[2] * sthe * sphi + Phi * dsthe[2] * sphi + Phi * sthe * dsphi[2];

    } else {
        CCTK_WARN(0, "Invalid choice for scalar field initial profile");
    }
}

/*-----------------------------------------------------------*/
