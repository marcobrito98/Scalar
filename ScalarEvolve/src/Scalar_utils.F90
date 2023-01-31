
#include "cctk.h"
#include "cctk_Parameters.h"

subroutine Scalar_apply_jacobian(dvar, jac)
  implicit none

  CCTK_REAL, intent(inout) :: dvar(3)
  CCTK_REAL, intent(in)    :: jac(3,3)
  CCTK_REAL                :: xdvar(3)
  CCTK_INT                 :: a, b

  xdvar = 0
  do a = 1, 3
     do b = 1, 3
        xdvar(a) = xdvar(a) + dvar(b) * jac(b,a)
     end do
  end do

  dvar = xdvar

end subroutine Scalar_apply_jacobian
!
!=============================================================================
!
subroutine Scalar_apply_jacobian2(dvar, ddvar, jac, hes)
  implicit none

  CCTK_REAL, intent(inout) :: ddvar(3,3), dvar(3)
  CCTK_REAL, intent(in)    :: jac(3,3), hes(3,3,3)
  CCTK_REAL                :: xddvar(3,3), xdvar(3)
  CCTK_INT                 :: a, b, c, d

  xdvar = 0
  do a = 1, 3
     do b = 1, 3
        xdvar(a) = xdvar(a) + dvar(b) * jac(b,a)
     end do
  end do

  xddvar = 0
  do a = 1, 3
     do b = 1, 3
        do c = 1, 3
           xddvar(a,b) = xddvar(a,b) + dvar(c) * hes(c,a,b)
           do d = 1, 3
              xddvar(a,b) = xddvar(a,b) + ddvar(c,d) * jac(c,a) * jac(d,b)
           end do
        end do
     end do
  end do

  dvar  = xdvar
  ddvar = xddvar

end subroutine Scalar_apply_jacobian2
