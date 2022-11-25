
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine Scalar_rhs_force( CCTK_ARGUMENTS )

 implicit none
 DECLARE_CCTK_ARGUMENTS
 DECLARE_CCTK_PARAMETERS

 CCTK_INT  i, j, k
 CCTK_REAL var1, var2                             ! Re and Im part of external forcing: to fill grid functions
 CCTK_REAL xp(3), tp                              ! grid point coordinates, time and local scalar field
 CCTK_REAL rr, rr2                                ! 3D radius
 CCTK_REAL rho, rho2                              ! 2D radius
 CCTK_REAL ctheta, ctheta2, cphi, cphi2, c2phi    ! Angles cos
 CCTK_REAL stheta, stheta2, sphi, sphi2, s2phi    ! Angles sin
 CCTK_REAL Pi 

 Pi = 4.0*ATAN(1.0)

  !$OMP PARALLEL DO COLLAPSE(3) &
  !$OMP PRIVATE(var1,var2,&
  !$OMP i, j, k,&
  !$OMP xp,tp,rr, rr2, rho, rho2,&
  !$OMP ctheta, ctheta2, cphi, cphi2, c2phi,&
  !$OMP stheta, stheta2, sphi, sphi2, s2phi)
  do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
  do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
  do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)

    xp(1) = x(i,j,k)
    xp(2) = y(i,j,k)
    xp(3) = z(i,j,k)

    tp = cctk_time

    rr2  = xp(1)*xp(1) + xp(2)*xp(2) + xp(3)*xp(3)
    rho2 = xp(1)*xp(1) + xp(2)*xp(2)

    if (rr2 < eps_r) then
        rr2 = eps_r
    end if

    if (rho2 < eps_r) then
        rho2 = eps_r
    end if

    rr   = SQRT(rr2)
    rho  = SQRT(rho2)

    ctheta  = xp(3) / rr;
    ctheta2 = ctheta * ctheta;
    cphi    = xp(1) / rho;
    cphi2   = cphi * cphi;

    stheta  = rho / rr;
    stheta2 = stheta * stheta;
    sphi    = xp(2) / rho;
    sphi2   = sphi * sphi;

    c2phi   = cphi2 - sphi2;
    s2phi   = 2.0 * cphi * sphi;


    if (lEF == 0 .and. mEF == 0) then
        var1 = 1.0 / SQRT( 4.0 * Pi )
        var2 = 0.0

    else if (lEF == 1 .and. mEF == 1 ) then
        var1 = - SQRT( 3.0 / ( 8.0 * Pi ) ) * cphi * stheta 
        var2 = - SQRT( 3.0 / ( 8.0 * Pi ) ) * sphi * stheta

    else if (lEF == 2 .and. mEF == 2 ) then
        var1 = SQRT( 15.0 / ( 32.0*Pi ) ) * c2phi * stheta2 
        var2 = SQRT( 15.0 / ( 32.0*Pi ) ) * s2phi * stheta2
    end if

    var1 = ampEF * EXP(- ( rr - rEF ) * ( rr - rEF ) / (widthEF * widthEF) ) * var1 * COS( omegaEF * tp )
    var2 = ampEF * EXP(- ( rr - rEF ) * ( rr - rEF ) / (widthEF * widthEF) ) * var2 * SIN( omegaEF * tp )

    !Fill grid function

    Fext1(i,j,k) = var1
    Fext2(i,j,k) = var2

       end do
     end do
   end do
   !$OMP END PARALLEL DO

end subroutine Scalar_rhs_force
