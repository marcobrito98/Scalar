
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine Scalar_ord4_calc_rhs( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  ! Fundamental variables
  CCTK_REAL                alph, beta(3)
  CCTK_REAL                ch, hh(3,3), hu(3,3), trk, dethh
  CCTK_REAL                lphi1, lphi2, lKphi1, lKphi2

  ! First derivatives
  CCTK_REAL                d1_alph(3)
  CCTK_REAL                d1_ch(3), d1_hh(3,3,3)
  CCTK_REAL                d1_lphi1(3), d1_lphi2(3), d1_lKphi1(3), d1_lKphi2(3)

  ! Second derivatives
  CCTK_REAL                d2_lphi1(3,3), d2_lphi2(3,3)

  ! Advection derivatives
  CCTK_REAL                ad1_lphi1, ad1_lphi2, ad1_lKphi1, ad1_lKphi2
  CCTK_REAL                d1_f(3)   ! Place holder for the advection derivs

  ! Covariant derivatives
  CCTK_REAL                cd2_lphi1(3,3), cd2_lphi2(3,3)

  ! Auxiliary variables
  CCTK_REAL                cf1(3,3,3), cf2(3,3,3)
  CCTK_REAL                tr_dalp_dphi1, tr_cd2_phi1, tr_dch_dphi1
  CCTK_REAL                tr_dalp_dphi2, tr_cd2_phi2, tr_dch_dphi2

  ! Right hand sides
  CCTK_REAL                rhs_lphi1, rhs_lphi2 
  CCTK_REAL                rhs_lKphi1, rhs_lKphi2

  ! Misc variables
  CCTK_REAL                dx12, dy12, dz12, dxsq12, dysq12, dzsq12,         &
                           dxdy144, dxdz144, dydz144
  CCTK_INT                 i, j, k
  CCTK_INT                 di, dj, dk
  CCTK_REAL, parameter ::  one = 1
  CCTK_INT                 a, b, c, m

  ! Excision variables
  CCTK_INT                 sn1, sn2
  CCTK_REAL                rsn1_2, rsn2_2, rr, lambda
  CCTK_REAL                Rout_excision1, Rout_excision2
  CCTK_REAL                Rin_excision1, Rin_excision2

  dx12 = 12*CCTK_DELTA_SPACE(1)
  dy12 = 12*CCTK_DELTA_SPACE(2)
  dz12 = 12*CCTK_DELTA_SPACE(3)

  dxsq12 = 12*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(1)
  dysq12 = 12*CCTK_DELTA_SPACE(2)*CCTK_DELTA_SPACE(2)
  dzsq12 = 12*CCTK_DELTA_SPACE(3)*CCTK_DELTA_SPACE(3)

  dxdy144 = 144*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2)
  dxdz144 = 144*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(3)
  dydz144 = 144*CCTK_DELTA_SPACE(2)*CCTK_DELTA_SPACE(3)

  rhs_phi1   = 0
  rhs_phi2   = 0
  rhs_Kphi1  = 0
  rhs_Kphi2  = 0

  ! convert ADM variables to BSSN-like ones
  call Scalar_adm2bssn(CCTK_PASS_FTOF)

  !$OMP PARALLEL DO COLLAPSE(3) &
  !$OMP PRIVATE(alph, beta, ch, hh, hu, trk, dethh,&
  !$OMP lphi1, lphi2, lKphi1, lKphi2,&
  !$OMP d1_alph, d1_ch, d1_hh,&
  !$OMP d1_lphi1, d1_lphi2, d1_lKphi1, d1_lKphi2,&
  !$OMP d2_lphi1, d2_lphi2, ad1_lphi1, ad1_lphi2, ad1_lKphi1, ad1_lKphi2,&
  !$OMP d1_f, cf1, cf2, cd2_lphi1, cd2_lphi2,&
  !$OMP tr_dalp_dphi1, tr_cd2_phi1, tr_dch_dphi1,&
  !$OMP tr_dalp_dphi2, tr_cd2_phi2, tr_dch_dphi2,&
  !$OMP rhs_lphi1, rhs_lphi2, rhs_lKphi1, rhs_lKphi2,&
  !$OMP i, j, k,&
  !$OMP di, dj, dk,&
  !$OMP a, b, c, m, rsn1_2, rsn2_2, rr, lambda)
  do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
  do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
  do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)

    !------------ Get local variables ----------

    ch      = chi(i,j,k)

    hh(1,1) = hxx(i,j,k)
    hh(1,2) = hxy(i,j,k)
    hh(1,3) = hxz(i,j,k)
    hh(2,2) = hyy(i,j,k)
    hh(2,3) = hyz(i,j,k)
    hh(3,3) = hzz(i,j,k)
    hh(2,1) = hh(1,2)
    hh(3,1) = hh(1,3)
    hh(3,2) = hh(2,3)

    trk     = tracek(i,j,k)

    alph    = alp(i,j,k)

    beta(1) = betax(i,j,k)
    beta(2) = betay(i,j,k)
    beta(3) = betaz(i,j,k)


    lphi1   = phi1(i,j,k)
    lphi2   = phi2(i,j,k)

    lKphi1  = Kphi1(i,j,k)
    lKphi2  = Kphi2(i,j,k)

    !-------------------------------------------


    !------------ Invert 3-metric ----------------
    ! NOTE: deth = 1 by construction, but that is not satisfied numerically
    dethh =       hh(1,1) * hh(2,2) * hh(3,3)                              &
            + 2 * hh(1,2) * hh(1,3) * hh(2,3)                              &
            -     hh(1,1) * hh(2,3) ** 2                                   &
            -     hh(2,2) * hh(1,3) ** 2                                   &
            -     hh(3,3) * hh(1,2) ** 2
    hu(1,1) = (hh(2,2) * hh(3,3) - hh(2,3) ** 2     ) / dethh
    hu(2,2) = (hh(1,1) * hh(3,3) - hh(1,3) ** 2     ) / dethh
    hu(3,3) = (hh(1,1) * hh(2,2) - hh(1,2) ** 2     ) / dethh
    hu(1,2) = (hh(1,3) * hh(2,3) - hh(1,2) * hh(3,3)) / dethh
    hu(1,3) = (hh(1,2) * hh(2,3) - hh(1,3) * hh(2,2)) / dethh
    hu(2,3) = (hh(1,3) * hh(1,2) - hh(2,3) * hh(1,1)) / dethh
    hu(2,1) = hu(1,2)
    hu(3,1) = hu(1,3)
    hu(3,2) = hu(2,3)
    !-------------------------------------------


    !------------ Centered 1st derivatives -----
    ! d1_ch(3)
    d1_ch(1) = (   -chi(i+2,j,k) + 8*chi(i+1,j,k)                          &
                - 8*chi(i-1,j,k) +   chi(i-2,j,k) ) / dx12

    d1_ch(2) = (   -chi(i,j+2,k) + 8*chi(i,j+1,k)                          &
                - 8*chi(i,j-1,k) +   chi(i,j-2,k) ) / dy12

    d1_ch(3) = (   -chi(i,j,k+2) + 8*chi(i,j,k+1)                          &
                - 8*chi(i,j,k-1) +   chi(i,j,k-2) ) / dz12

    ! d1_hh(3,3,3)
    d1_hh(1,1,1) = (   -hxx(i+2,j,k) + 8*hxx(i+1,j,k)                      &
                    - 8*hxx(i-1,j,k) +   hxx(i-2,j,k) ) / dx12
    d1_hh(1,2,1) = (   -hxy(i+2,j,k) + 8*hxy(i+1,j,k)                      &
                    - 8*hxy(i-1,j,k) +   hxy(i-2,j,k) ) / dx12
    d1_hh(1,3,1) = (   -hxz(i+2,j,k) + 8*hxz(i+1,j,k)                      &
                    - 8*hxz(i-1,j,k) +   hxz(i-2,j,k) ) / dx12
    d1_hh(2,2,1) = (   -hyy(i+2,j,k) + 8*hyy(i+1,j,k)                      &
                    - 8*hyy(i-1,j,k) +   hyy(i-2,j,k) ) / dx12
    d1_hh(2,3,1) = (   -hyz(i+2,j,k) + 8*hyz(i+1,j,k)                      &
                    - 8*hyz(i-1,j,k) +   hyz(i-2,j,k) ) / dx12
    d1_hh(3,3,1) = (   -hzz(i+2,j,k) + 8*hzz(i+1,j,k)                      &
                    - 8*hzz(i-1,j,k) +   hzz(i-2,j,k) ) / dx12

    d1_hh(1,1,2) = (   -hxx(i,j+2,k) + 8*hxx(i,j+1,k)                      &
                    - 8*hxx(i,j-1,k) +   hxx(i,j-2,k) ) / dy12
    d1_hh(1,2,2) = (   -hxy(i,j+2,k) + 8*hxy(i,j+1,k)                      &
                    - 8*hxy(i,j-1,k) +   hxy(i,j-2,k) ) / dy12
    d1_hh(1,3,2) = (   -hxz(i,j+2,k) + 8*hxz(i,j+1,k)                      &
                    - 8*hxz(i,j-1,k) +   hxz(i,j-2,k) ) / dy12
    d1_hh(2,2,2) = (   -hyy(i,j+2,k) + 8*hyy(i,j+1,k)                      &
                    - 8*hyy(i,j-1,k) +   hyy(i,j-2,k) ) / dy12
    d1_hh(2,3,2) = (   -hyz(i,j+2,k) + 8*hyz(i,j+1,k)                      &
                    - 8*hyz(i,j-1,k) +   hyz(i,j-2,k) ) / dy12
    d1_hh(3,3,2) = (   -hzz(i,j+2,k) + 8*hzz(i,j+1,k)                      &
                    - 8*hzz(i,j-1,k) +   hzz(i,j-2,k) ) / dy12

    d1_hh(1,1,3) = (   -hxx(i,j,k+2) + 8*hxx(i,j,k+1)                      &
                    - 8*hxx(i,j,k-1) +   hxx(i,j,k-2) ) / dz12
    d1_hh(1,2,3) = (   -hxy(i,j,k+2) + 8*hxy(i,j,k+1)                      &
                    - 8*hxy(i,j,k-1) +   hxy(i,j,k-2) ) / dz12
    d1_hh(1,3,3) = (   -hxz(i,j,k+2) + 8*hxz(i,j,k+1)                      &
                    - 8*hxz(i,j,k-1) +   hxz(i,j,k-2) ) / dz12
    d1_hh(2,2,3) = (   -hyy(i,j,k+2) + 8*hyy(i,j,k+1)                      &
                    - 8*hyy(i,j,k-1) +   hyy(i,j,k-2) ) / dz12
    d1_hh(2,3,3) = (   -hyz(i,j,k+2) + 8*hyz(i,j,k+1)                      &
                    - 8*hyz(i,j,k-1) +   hyz(i,j,k-2) ) / dz12
    d1_hh(3,3,3) = (   -hzz(i,j,k+2) + 8*hzz(i,j,k+1)                      &
                    - 8*hzz(i,j,k-1) +   hzz(i,j,k-2) ) / dz12

    d1_hh(2,1,:) = d1_hh(1,2,:)
    d1_hh(3,1,:) = d1_hh(1,3,:)
    d1_hh(3,2,:) = d1_hh(2,3,:)


    ! d1_alph(3)
    d1_alph(1) = (   -alp(i+2,j,k) + 8*alp(i+1,j,k)                        &
                  - 8*alp(i-1,j,k) +   alp(i-2,j,k) ) / dx12

    d1_alph(2) = (   -alp(i,j+2,k) + 8*alp(i,j+1,k)                        &
                  - 8*alp(i,j-1,k) +   alp(i,j-2,k) ) / dy12

    d1_alph(3) = (   -alp(i,j,k+2) + 8*alp(i,j,k+1)                        &
                  - 8*alp(i,j,k-1) +   alp(i,j,k-2) ) / dz12


    ! d1_lphi1(3)
    d1_lphi1(1) = (   -phi1(i+2,j,k) + 8*phi1(i+1,j,k)                          &
                   - 8*phi1(i-1,j,k) +   phi1(i-2,j,k) ) / dx12

    d1_lphi1(2) = (   -phi1(i,j+2,k) + 8*phi1(i,j+1,k)                          &
                   - 8*phi1(i,j-1,k) +   phi1(i,j-2,k) ) / dy12

    d1_lphi1(3) = (   -phi1(i,j,k+2) + 8*phi1(i,j,k+1)                          &
                   - 8*phi1(i,j,k-1) +   phi1(i,j,k-2) ) / dz12

    ! d1_lphi2(3)
    d1_lphi2(1) = (   -phi2(i+2,j,k) + 8*phi2(i+1,j,k)                          &
                   - 8*phi2(i-1,j,k) +   phi2(i-2,j,k) ) / dx12

    d1_lphi2(2) = (   -phi2(i,j+2,k) + 8*phi2(i,j+1,k)                          &
                   - 8*phi2(i,j-1,k) +   phi2(i,j-2,k) ) / dy12

    d1_lphi2(3) = (   -phi2(i,j,k+2) + 8*phi2(i,j,k+1)                          &
                   - 8*phi2(i,j,k-1) +   phi2(i,j,k-2) ) / dz12

    ! d1_lKphi1(3)
    d1_lKphi1(1) = (   -Kphi1(i+2,j,k) + 8*Kphi1(i+1,j,k)                          &
                    - 8*Kphi1(i-1,j,k) +   Kphi1(i-2,j,k) ) / dx12

    d1_lKphi1(2) = (   -Kphi1(i,j+2,k) + 8*Kphi1(i,j+1,k)                          &
                    - 8*Kphi1(i,j-1,k) +   Kphi1(i,j-2,k) ) / dy12

    d1_lKphi1(3) = (   -Kphi1(i,j,k+2) + 8*Kphi1(i,j,k+1)                          &
                    - 8*Kphi1(i,j,k-1) +   Kphi1(i,j,k-2) ) / dz12

    ! d1_lKphi2(3)
    d1_lKphi2(1) = (   -Kphi2(i+2,j,k) + 8*Kphi2(i+1,j,k)                          &
                    - 8*Kphi2(i-1,j,k) +   Kphi2(i-2,j,k) ) / dx12

    d1_lKphi2(2) = (   -Kphi2(i,j+2,k) + 8*Kphi2(i,j+1,k)                          &
                    - 8*Kphi2(i,j-1,k) +   Kphi2(i,j-2,k) ) / dy12

    d1_lKphi2(3) = (   -Kphi2(i,j,k+2) + 8*Kphi2(i,j,k+1)                          &
                    - 8*Kphi2(i,j,k-1) +   Kphi2(i,j,k-2) ) / dz12


    !--------------------------------------------------




    !------------- Centered 2nd derivatives -----------

    ! d2_lphi1(3,3)
    d2_lphi1(1,1) = (   -phi1(i+2,j,k) + 16*phi1(i+1,j,k) - 30*phi1(i,j,k)     &
                    + 16*phi1(i-1,j,k) -    phi1(i-2,j,k) ) / dxsq12

    d2_lphi1(2,2) = (   -phi1(i,j+2,k) + 16*phi1(i,j+1,k) - 30*phi1(i,j,k)     &
                    + 16*phi1(i,j-1,k) -    phi1(i,j-2,k) ) / dysq12

    d2_lphi1(3,3) = (   -phi1(i,j,k+2) + 16*phi1(i,j,k+1) - 30*phi1(i,j,k)     &
                    + 16*phi1(i,j,k-1) -    phi1(i,j,k-2) ) / dzsq12

    d2_lphi1(1,2) = (   -phi1(i-2,j+2,k) +  8*phi1(i-1,j+2,k) -  8*phi1(i+1,j+2,k) +   phi1(i+2,j+2,k)   &
                     + 8*phi1(i-2,j+1,k) - 64*phi1(i-1,j+1,k) + 64*phi1(i+1,j+1,k) - 8*phi1(i+2,j+1,k)   &
                     - 8*phi1(i-2,j-1,k) + 64*phi1(i-1,j-1,k) - 64*phi1(i+1,j-1,k) + 8*phi1(i+2,j-1,k)   &
                     +   phi1(i-2,j-2,k) -  8*phi1(i-1,j-2,k) +  8*phi1(i+1,j-2,k) -   phi1(i+2,j-2,k) ) / dxdy144

    d2_lphi1(1,3) = (   -phi1(i-2,j,k+2) +  8*phi1(i-1,j,k+2) -  8*phi1(i+1,j,k+2) +   phi1(i+2,j,k+2)   &
                     + 8*phi1(i-2,j,k+1) - 64*phi1(i-1,j,k+1) + 64*phi1(i+1,j,k+1) - 8*phi1(i+2,j,k+1)   &
                     - 8*phi1(i-2,j,k-1) + 64*phi1(i-1,j,k-1) - 64*phi1(i+1,j,k-1) + 8*phi1(i+2,j,k-1)   &
                     +   phi1(i-2,j,k-2) -  8*phi1(i-1,j,k-2) +  8*phi1(i+1,j,k-2) -   phi1(i+2,j,k-2) ) / dxdz144

    d2_lphi1(2,3) = (   -phi1(i,j-2,k+2) +  8*phi1(i,j-1,k+2) -  8*phi1(i,j+1,k+2) +   phi1(i,j+2,k+2)   &
                     + 8*phi1(i,j-2,k+1) - 64*phi1(i,j-1,k+1) + 64*phi1(i,j+1,k+1) - 8*phi1(i,j+2,k+1)   &
                     - 8*phi1(i,j-2,k-1) + 64*phi1(i,j-1,k-1) - 64*phi1(i,j+1,k-1) + 8*phi1(i,j+2,k-1)   &
                     +   phi1(i,j-2,k-2) -  8*phi1(i,j-1,k-2) +  8*phi1(i,j+1,k-2) -   phi1(i,j+2,k-2) ) / dydz144

    d2_lphi1(2,1) = d2_lphi1(1,2)
    d2_lphi1(3,1) = d2_lphi1(1,3)
    d2_lphi1(3,2) = d2_lphi1(2,3)

    ! d2_lphi2(3,3)
    d2_lphi2(1,1) = (   -phi2(i+2,j,k) + 16*phi2(i+1,j,k) - 30*phi2(i,j,k)     &
                    + 16*phi2(i-1,j,k) -    phi2(i-2,j,k) ) / dxsq12

    d2_lphi2(2,2) = (   -phi2(i,j+2,k) + 16*phi2(i,j+1,k) - 30*phi2(i,j,k)     &
                    + 16*phi2(i,j-1,k) -    phi2(i,j-2,k) ) / dysq12

    d2_lphi2(3,3) = (   -phi2(i,j,k+2) + 16*phi2(i,j,k+1) - 30*phi2(i,j,k)     &
                    + 16*phi2(i,j,k-1) -    phi2(i,j,k-2) ) / dzsq12

    d2_lphi2(1,2) = (   -phi2(i-2,j+2,k) +  8*phi2(i-1,j+2,k) -  8*phi2(i+1,j+2,k) +   phi2(i+2,j+2,k)   &
                     + 8*phi2(i-2,j+1,k) - 64*phi2(i-1,j+1,k) + 64*phi2(i+1,j+1,k) - 8*phi2(i+2,j+1,k)   &
                     - 8*phi2(i-2,j-1,k) + 64*phi2(i-1,j-1,k) - 64*phi2(i+1,j-1,k) + 8*phi2(i+2,j-1,k)   &
                     +   phi2(i-2,j-2,k) -  8*phi2(i-1,j-2,k) +  8*phi2(i+1,j-2,k) -   phi2(i+2,j-2,k) ) / dxdy144

    d2_lphi2(1,3) = (   -phi2(i-2,j,k+2) +  8*phi2(i-1,j,k+2) -  8*phi2(i+1,j,k+2) +   phi2(i+2,j,k+2)   &
                     + 8*phi2(i-2,j,k+1) - 64*phi2(i-1,j,k+1) + 64*phi2(i+1,j,k+1) - 8*phi2(i+2,j,k+1)   &
                     - 8*phi2(i-2,j,k-1) + 64*phi2(i-1,j,k-1) - 64*phi2(i+1,j,k-1) + 8*phi2(i+2,j,k-1)   &
                     +   phi2(i-2,j,k-2) -  8*phi2(i-1,j,k-2) +  8*phi2(i+1,j,k-2) -   phi2(i+2,j,k-2) ) / dxdz144

    d2_lphi2(2,3) = (   -phi2(i,j-2,k+2) +  8*phi2(i,j-1,k+2) -  8*phi2(i,j+1,k+2) +   phi2(i,j+2,k+2)   &
                     + 8*phi2(i,j-2,k+1) - 64*phi2(i,j-1,k+1) + 64*phi2(i,j+1,k+1) - 8*phi2(i,j+2,k+1)   &
                     - 8*phi2(i,j-2,k-1) + 64*phi2(i,j-1,k-1) - 64*phi2(i,j+1,k-1) + 8*phi2(i,j+2,k-1)   &
                     +   phi2(i,j-2,k-2) -  8*phi2(i,j-1,k-2) +  8*phi2(i,j+1,k-2) -   phi2(i,j+2,k-2) ) / dydz144

    d2_lphi2(2,1) = d2_lphi2(1,2)
    d2_lphi2(3,1) = d2_lphi2(1,3)
    d2_lphi2(3,2) = d2_lphi2(2,3)


    !------------ Advection derivatives --------
    if( use_advection_stencils /= 0 ) then

      di = int( sign( one, beta(1) ) )
      dj = int( sign( one, beta(2) ) )
      dk = int( sign( one, beta(3) ) )

      ! ad1_lphi1
      d1_f(1) = di * ( -3*phi1(i-di,j,k) - 10*phi1(i,j,k) + 18*phi1(i+di,j,k)   &
                      - 6*phi1(i+2*di,j,k)  + phi1(i+3*di,j,k)) / dx12
      d1_f(2) = dj * ( -3*phi1(i,j-dj,k) - 10*phi1(i,j,k) + 18*phi1(i,j+dj,k)   &
                      - 6*phi1(i,j+2*dj,k)  + phi1(i,j+3*dj,k)) / dy12
      d1_f(3) = dk * ( -3*phi1(i,j,k-dk) - 10*phi1(i,j,k) + 18*phi1(i,j,k+dk)   &
                      - 6*phi1(i,j,k+2*dk)  + phi1(i,j,k+3*dk)) / dz12
      ad1_lphi1 = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

      ! ad1_lKphi1
      d1_f(1) = di * ( -3*Kphi1(i-di,j,k) - 10*Kphi1(i,j,k) + 18*Kphi1(i+di,j,k)   &
                      - 6*Kphi1(i+2*di,j,k)  + Kphi1(i+3*di,j,k)) / dx12
      d1_f(2) = dj * ( -3*Kphi1(i,j-dj,k) - 10*Kphi1(i,j,k) + 18*Kphi1(i,j+dj,k)   &
                      - 6*Kphi1(i,j+2*dj,k)  + Kphi1(i,j+3*dj,k)) / dy12
      d1_f(3) = dk * ( -3*Kphi1(i,j,k-dk) - 10*Kphi1(i,j,k) + 18*Kphi1(i,j,k+dk)   &
                      - 6*Kphi1(i,j,k+2*dk)  + Kphi1(i,j,k+3*dk)) / dz12
      ad1_lKphi1 = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

      ! ad1_lphi2
      d1_f(1) = di * ( -3*phi2(i-di,j,k) - 10*phi2(i,j,k) + 18*phi2(i+di,j,k)   &
                      - 6*phi2(i+2*di,j,k)  + phi2(i+3*di,j,k)) / dx12
      d1_f(2) = dj * ( -3*phi2(i,j-dj,k) - 10*phi2(i,j,k) + 18*phi2(i,j+dj,k)   &
                      - 6*phi2(i,j+2*dj,k)  + phi2(i,j+3*dj,k)) / dy12
      d1_f(3) = dk * ( -3*phi2(i,j,k-dk) - 10*phi2(i,j,k) + 18*phi2(i,j,k+dk)   &
                      - 6*phi2(i,j,k+2*dk)  + phi2(i,j,k+3*dk)) / dz12
      ad1_lphi2 = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

      ! ad1_lKphi2
      d1_f(1) = di * ( -3*Kphi2(i-di,j,k) - 10*Kphi2(i,j,k) + 18*Kphi2(i+di,j,k)   &
                      - 6*Kphi2(i+2*di,j,k)  + Kphi2(i+3*di,j,k)) / dx12
      d1_f(2) = dj * ( -3*Kphi2(i,j-dj,k) - 10*Kphi2(i,j,k) + 18*Kphi2(i,j+dj,k)   &
                      - 6*Kphi2(i,j+2*dj,k)  + Kphi2(i,j+3*dj,k)) / dy12
      d1_f(3) = dk * ( -3*Kphi2(i,j,k-dk) - 10*Kphi2(i,j,k) + 18*Kphi2(i,j,k+dk)   &
                      - 6*Kphi2(i,j,k+2*dk)  + Kphi2(i,j,k+3*dk)) / dz12
      ad1_lKphi2 = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

    else

      ! ad1_lphi1
      ad1_lphi1 = beta(1)*d1_lphi1(1) + beta(2)*d1_lphi1(2) + beta(3)*d1_lphi1(3)

      ! ad1_lKphi1
      ad1_lKphi1 = beta(1)*d1_lKphi1(1) + beta(2)*d1_lKphi1(2) + beta(3)*d1_lKphi1(3)

      ! ad1_lphi2
      ad1_lphi2 = beta(1)*d1_lphi2(1) + beta(2)*d1_lphi2(2) + beta(3)*d1_lphi2(3)

      ! ad1_lKphi2
      ad1_lKphi2 = beta(1)*d1_lKphi2(1) + beta(2)*d1_lKphi2(2) + beta(3)*d1_lKphi2(3)

    end if

    !-------------------------------------------


    !------------ Christoffel symbols ----------
    cf1 = 0
    do a = 1, 3
      do b = 1, 3
        do c = b, 3
          cf1(a,b,c) = 0.5d0 * (d1_hh(a,b,c) + d1_hh(a,c,b) - d1_hh(b,c,a))
        end do
      end do
    end do
    cf1(:,2,1) = cf1(:,1,2)
    cf1(:,3,1) = cf1(:,1,3)
    cf1(:,3,2) = cf1(:,2,3)

    cf2 = 0
    do a = 1, 3
      do b = 1, 3
        do c = b, 3
          do m = 1, 3
            cf2(a,b,c) = cf2(a,b,c) + hu(a,m) * cf1(m,b,c)
          end do
        end do
      end do
    end do
    cf2(:,2,1) = cf2(:,1,2)
    cf2(:,3,1) = cf2(:,1,3)
    cf2(:,3,2) = cf2(:,2,3)
    !-------------------------------------------


    !------------ Covariant derivatives --------
    cd2_lphi1  = d2_lphi1
    cd2_lphi2  = d2_lphi2
    do a = 1, 3
      do b = a, 3
        do m = 1, 3
          cd2_lphi1(a,b)  = cd2_lphi1(a,b) - cf2(m,a,b) * d1_lphi1(m)
          cd2_lphi2(a,b)  = cd2_lphi2(a,b) - cf2(m,a,b) * d1_lphi2(m)
        end do
      end do
    end do
    cd2_lphi1(2,1) = cd2_lphi1(1,2)
    cd2_lphi1(3,1) = cd2_lphi1(1,3)
    cd2_lphi1(3,2) = cd2_lphi1(2,3)
    cd2_lphi2(2,1) = cd2_lphi2(1,2)
    cd2_lphi2(3,1) = cd2_lphi2(1,3)
    cd2_lphi2(3,2) = cd2_lphi2(2,3)

    !-------------------------------------------

    !------------ Advection and Twist terms ----

    ! rhs_lphi1, rhs_lphi2
    rhs_lphi1  = ad1_lphi1
    rhs_lphi2  = ad1_lphi2

    ! rhs_lKphi1, rhs_lKphi2
    rhs_lKphi1 = ad1_lKphi1
    rhs_lKphi2 = ad1_lKphi2

    !-------------------------------------------


    !------------ Source terms -----------------

    ! rhs_lphi1, rhs_lphi2 
    rhs_lphi1 = rhs_lphi1 - 2.0 * alph * lKphi1

    rhs_lphi2 = rhs_lphi2 - 2.0 * alph * lKphi2


    ! rhs_lKphi1, rhs_lKphi2
    tr_dalp_dphi1 = 0
    tr_dalp_dphi2 = 0
    tr_cd2_phi1   = 0
    tr_cd2_phi2   = 0
    tr_dch_dphi1  = 0
    tr_dch_dphi2  = 0
    do a = 1, 3
      do b = 1, 3
        tr_dalp_dphi1 = tr_dalp_dphi1 + hu(a,b) * d1_alph(a) * d1_lphi1(b)
        tr_dalp_dphi2 = tr_dalp_dphi2 + hu(a,b) * d1_alph(a) * d1_lphi2(b)
        tr_cd2_phi1   = tr_cd2_phi1   + hu(a,b) * cd2_lphi1(a,b)
        tr_cd2_phi2   = tr_cd2_phi2   + hu(a,b) * cd2_lphi2(a,b)
        tr_dch_dphi1  = tr_dch_dphi1  + hu(a,b) * d1_ch(a) * d1_lphi1(b)
        tr_dch_dphi2  = tr_dch_dphi2  + hu(a,b) * d1_ch(a) * d1_lphi2(b)
      end do
    end do

    rhs_lKphi1 = rhs_lKphi1 - 0.5 * ch * tr_dalp_dphi1                         &
                 + 0.5 * alph * ( - ch * tr_cd2_phi1 + 0.5 * tr_dch_dphi1      &
                                  + mu*mu * lphi1 + 2.0 * trk * lKphi1 )

    rhs_lKphi2 = rhs_lKphi2 - 0.5 * ch * tr_dalp_dphi2                         &
                 + 0.5 * alph * ( - ch * tr_cd2_phi2 + 0.5 * tr_dch_dphi2      &
                                  + mu*mu * lphi2 + 2.0 * trk * lKphi2 )

    !-------------------------------------------

    !------------ Excise and write to grid functions ------
    if (use_excision /= 0) then
       sn1 = excision_surface(1) + 1
       sn2 = excision_surface(2) + 1

       ! write(*,*) "sf_active = ", sf_active(sn1), sf_active(sn2)
       ! write(*,*) "sf_mean_radius = ", sf_mean_radius(sn1), sf_mean_radius(sn2)
       ! write(*,*) "sf_origin (sn1) = ", sf_origin_x(sn1), sf_origin_y(sn1), sf_origin_z(sn1)
       ! write(*,*) "sf_origin (sn2) = ", sf_origin_x(sn2), sf_origin_y(sn2), sf_origin_z(sn2)

       if (sf_active(sn1) /= 0 .and. sf_active(sn2) /= 0) then
          Rout_excision1 = excision_surface_rout_factor(1) * sf_mean_radius(sn1)
          Rout_excision2 = excision_surface_rout_factor(2) * sf_mean_radius(sn2)
          Rin_excision1  = excision_surface_rin_factor(1)  * sf_mean_radius(sn1)
          Rin_excision2  = excision_surface_rin_factor(2)  * sf_mean_radius(sn2)

          rsn1_2 =  (x(i,j,k) - sf_origin_x(sn1))**2      &
                  + (y(i,j,k) - sf_origin_y(sn1))**2      &
                  + (z(i,j,k) - sf_origin_z(sn1))**2

          rsn2_2 =  (x(i,j,k) - sf_origin_x(sn2))**2      &
                  + (y(i,j,k) - sf_origin_y(sn2))**2      &
                  + (z(i,j,k) - sf_origin_z(sn2))**2
       else
          call CCTK_ERROR("Using excision, but surface is not active!")
       end if

       if (rsn1_2 > Rout_excision1 * Rout_excision1 .and.       &
           rsn2_2 > Rout_excision2 * Rout_excision2) then
          ! regular points, outside horizon

          rhs_phi1(i,j,k)  = rhs_lphi1
          rhs_phi2(i,j,k)  = rhs_lphi2
          rhs_Kphi1(i,j,k) = rhs_lKphi1
          rhs_Kphi2(i,j,k) = rhs_lKphi2
       else if(rsn1_2 < Rin_excision1 * Rin_excision1 .or.      &
               rsn2_2 < Rin_excision2 * Rin_excision2) then
          ! excision zone
          rhs_phi1(i,j,k)  = 0.0
          rhs_phi2(i,j,k)  = 0.0
          rhs_Kphi1(i,j,k) = 0.0
          rhs_Kphi2(i,j,k) = 0.0
       else
          ! buffer zone
          if (rsn1_2 < rsn2_2) then
             rr = sqrt(rsn1_2)
             lambda = 0.5 * (1 + tanh( 1./(Rin_excision1 - rr) - 1./(rr - Rout_excision1) ) )
          else
             rr = sqrt(rsn2_2)
             lambda = 0.5 * (1 + tanh( 1./(Rin_excision2 - rr) - 1./(rr - Rout_excision2) ) )
          end if
          rhs_phi1(i,j,k)  = lambda * rhs_lphi1
          rhs_phi2(i,j,k)  = lambda * rhs_lphi2
          rhs_Kphi1(i,j,k) = lambda * rhs_lKphi1
          rhs_Kphi2(i,j,k) = lambda * rhs_lKphi2
       end if

    else ! no excision
       rhs_phi1(i,j,k)  = rhs_lphi1
       rhs_phi2(i,j,k)  = rhs_lphi2
       rhs_Kphi1(i,j,k) = rhs_lKphi1
       rhs_Kphi2(i,j,k) = rhs_lKphi2
    end if
    !-------------------------------------------

    !if(i == 51 .and. j == 5 .and. k == 15) then
    !  !if(cur_it == 2) then
    !  write(*,*) 'cur_it = ', cur_it
    !  write(*,*) x(i,j,k), y(i,j,k), z(i,j,k)
    !  !end if
    !  write(*,*) ''
    !  write(*,*)
    !  call flush(6)
    !  !if( cur_it == 2 ) call CCTK_WARN( 0, '==============' )
    !end if

  end do
  end do
  end do
  !$OMP END PARALLEL DO

end subroutine Scalar_ord4_calc_rhs
!
!=============================================================================
!
subroutine Scalar_ord4_calc_rhs_bdry( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL, parameter :: one = 1.0
  CCTK_INT ierr

  ierr = NewRad_Apply(cctkGH, phi1, rhs_phi1, phi1_0, one, n_phi1)
  ierr = NewRad_Apply(cctkGH, phi2, rhs_phi2, phi2_0, one, n_phi2)

  ierr = NewRad_Apply(cctkGH, Kphi1, rhs_Kphi1, Kphi1_0, one, n_Kphi1)
  ierr = NewRad_Apply(cctkGH, Kphi2, rhs_Kphi2, Kphi2_0, one, n_Kphi2)

end subroutine Scalar_ord4_calc_rhs_bdry
