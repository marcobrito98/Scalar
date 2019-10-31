
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine Scalar_ord4_calc_rhs( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  ! Fundamental variables
  CCTK_REAL                alph, beta(3), beta_l(3)    ! beta is beta^{i} (upper index)
  CCTK_REAL                ch, hh(3,3), hu(3,3), trk, dethh
  CCTK_REAL                lphi1, lphi2, lKphi1, lKphi2

  ! First derivatives
  CCTK_REAL                d1_hh11(3), d1_hh12(3), d1_hh13(3), d1_hh22(3), d1_hh23(3), d1_hh33(3)
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

  ! jacobian
  integer                  istat
  logical                  use_jacobian
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ11, lJ12, lJ13
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ21, lJ22, lJ23
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ31, lJ32, lJ33
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: ldJ111, ldJ112, ldJ113, ldJ122, ldJ123, ldJ133
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: ldJ211, ldJ212, ldJ213, ldJ222, ldJ223, ldJ233
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: ldJ311, ldJ312, ldJ313, ldJ322, ldJ323, ldJ333

  CCTK_POINTER             lJ11_ptr, lJ12_ptr, lJ13_ptr
  CCTK_POINTER             lJ21_ptr, lJ22_ptr, lJ23_ptr
  CCTK_POINTER             lJ31_ptr, lJ32_ptr, lJ33_ptr
  CCTK_POINTER             ldJ111_ptr, ldJ112_ptr, ldJ113_ptr, ldJ122_ptr, ldJ123_ptr, ldJ133_ptr
  CCTK_POINTER             ldJ211_ptr, ldJ212_ptr, ldJ213_ptr, ldJ222_ptr, ldJ223_ptr, ldJ233_ptr
  CCTK_POINTER             ldJ311_ptr, ldJ312_ptr, ldJ313_ptr, ldJ322_ptr, ldJ323_ptr, ldJ333_ptr

  CCTK_REAL                jac(3,3), hes(3,3,3)

  pointer (lJ11_ptr, lJ11), (lJ12_ptr, lJ12), (lJ13_ptr, lJ13)
  pointer (lJ21_ptr, lJ21), (lJ22_ptr, lJ22), (lJ23_ptr, lJ23)
  pointer (lJ31_ptr, lJ31), (lJ32_ptr, lJ32), (lJ33_ptr, lJ33)

  pointer (ldJ111_ptr, ldJ111), (ldJ112_ptr, ldJ112), (ldJ113_ptr, ldJ113), (ldJ122_ptr, ldJ122), (ldJ123_ptr, ldJ123), (ldJ133_ptr, ldJ133)
  pointer (ldJ211_ptr, ldJ211), (ldJ212_ptr, ldJ212), (ldJ213_ptr, ldJ213), (ldJ222_ptr, ldJ222), (ldJ223_ptr, ldJ223), (ldJ233_ptr, ldJ233)
  pointer (ldJ311_ptr, ldJ311), (ldJ312_ptr, ldJ312), (ldJ313_ptr, ldJ313), (ldJ322_ptr, ldJ322), (ldJ323_ptr, ldJ323), (ldJ333_ptr, ldJ333)

  ! TODO: can this be active but with a cartesian mapping choice?
  call CCTK_IsFunctionAliased(istat, "MultiPatch_GetDomainSpecification")
  if (istat == 0) then
     use_jacobian = .false.
  else
     use_jacobian = .true.
  end if

  if (use_jacobian) then
     call CCTK_VarDataPtr(lJ11_ptr, cctkGH, 0, "Coordinates::J11")
     call CCTK_VarDataPtr(lJ12_ptr, cctkGH, 0, "Coordinates::J12")
     call CCTK_VarDataPtr(lJ13_ptr, cctkGH, 0, "Coordinates::J13")
     call CCTK_VarDataPtr(lJ21_ptr, cctkGH, 0, "Coordinates::J21")
     call CCTK_VarDataPtr(lJ22_ptr, cctkGH, 0, "Coordinates::J22")
     call CCTK_VarDataPtr(lJ23_ptr, cctkGH, 0, "Coordinates::J23")
     call CCTK_VarDataPtr(lJ31_ptr, cctkGH, 0, "Coordinates::J31")
     call CCTK_VarDataPtr(lJ32_ptr, cctkGH, 0, "Coordinates::J32")
     call CCTK_VarDataPtr(lJ33_ptr, cctkGH, 0, "Coordinates::J33")

     call CCTK_VarDataPtr(ldJ111_ptr, cctkGH, 0, "Coordinates::dJ111")
     call CCTK_VarDataPtr(ldJ112_ptr, cctkGH, 0, "Coordinates::dJ112")
     call CCTK_VarDataPtr(ldJ113_ptr, cctkGH, 0, "Coordinates::dJ113")
     call CCTK_VarDataPtr(ldJ122_ptr, cctkGH, 0, "Coordinates::dJ122")
     call CCTK_VarDataPtr(ldJ123_ptr, cctkGH, 0, "Coordinates::dJ123")
     call CCTK_VarDataPtr(ldJ133_ptr, cctkGH, 0, "Coordinates::dJ133")

     call CCTK_VarDataPtr(ldJ211_ptr, cctkGH, 0, "Coordinates::dJ211")
     call CCTK_VarDataPtr(ldJ212_ptr, cctkGH, 0, "Coordinates::dJ212")
     call CCTK_VarDataPtr(ldJ213_ptr, cctkGH, 0, "Coordinates::dJ213")
     call CCTK_VarDataPtr(ldJ222_ptr, cctkGH, 0, "Coordinates::dJ222")
     call CCTK_VarDataPtr(ldJ223_ptr, cctkGH, 0, "Coordinates::dJ223")
     call CCTK_VarDataPtr(ldJ233_ptr, cctkGH, 0, "Coordinates::dJ233")

     call CCTK_VarDataPtr(ldJ311_ptr, cctkGH, 0, "Coordinates::dJ311")
     call CCTK_VarDataPtr(ldJ312_ptr, cctkGH, 0, "Coordinates::dJ312")
     call CCTK_VarDataPtr(ldJ313_ptr, cctkGH, 0, "Coordinates::dJ313")
     call CCTK_VarDataPtr(ldJ322_ptr, cctkGH, 0, "Coordinates::dJ322")
     call CCTK_VarDataPtr(ldJ323_ptr, cctkGH, 0, "Coordinates::dJ323")
     call CCTK_VarDataPtr(ldJ333_ptr, cctkGH, 0, "Coordinates::dJ333")
end if


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
  !$OMP d1_hh11, d1_hh12, d1_hh13, d1_hh22, d1_hh23, d1_hh33,&
  !$OMP d1_alph, d1_ch, d1_hh,&
  !$OMP d1_lphi1, d1_lphi2, d1_lKphi1, d1_lKphi2,&
  !$OMP d2_lphi1, d2_lphi2, ad1_lphi1, ad1_lphi2, ad1_lKphi1, ad1_lKphi2,&
  !$OMP d1_f, cf1, cf2, cd2_lphi1, cd2_lphi2,&
  !$OMP tr_dalp_dphi1, tr_cd2_phi1, tr_dch_dphi1,&
  !$OMP tr_dalp_dphi2, tr_cd2_phi2, tr_dch_dphi2,&
  !$OMP rhs_lphi1, rhs_lphi2, rhs_lKphi1, rhs_lKphi2,&
  !$OMP i, j, k,&
  !$OMP di, dj, dk,&
  !$OMP a, b, c, m, rsn1_2, rsn2_2, rr, lambda, &
  !$OMP jac, hes, beta_l)
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

    ! hes(i,j,k) = jac(i,j),k
    if (use_jacobian) then
       jac(1,1) = lJ11(i,j,k)
       jac(1,2) = lJ12(i,j,k)
       jac(1,3) = lJ13(i,j,k)
       jac(2,1) = lJ21(i,j,k)
       jac(2,2) = lJ22(i,j,k)
       jac(2,3) = lJ23(i,j,k)
       jac(3,1) = lJ31(i,j,k)
       jac(3,2) = lJ32(i,j,k)
       jac(3,3) = lJ33(i,j,k)

       hes(1,1,1) = ldJ111(i,j,k)
       hes(1,1,2) = ldJ112(i,j,k)
       hes(1,1,3) = ldJ113(i,j,k)
       hes(1,2,1) = ldJ112(i,j,k)
       hes(1,2,2) = ldJ122(i,j,k)
       hes(1,2,3) = ldJ123(i,j,k)
       hes(1,3,1) = ldJ113(i,j,k)
       hes(1,3,2) = ldJ123(i,j,k)
       hes(1,3,3) = ldJ133(i,j,k)

       hes(2,1,1) = ldJ211(i,j,k)
       hes(2,1,2) = ldJ212(i,j,k)
       hes(2,1,3) = ldJ213(i,j,k)
       hes(2,2,1) = ldJ212(i,j,k)
       hes(2,2,2) = ldJ222(i,j,k)
       hes(2,2,3) = ldJ223(i,j,k)
       hes(2,3,1) = ldJ213(i,j,k)
       hes(2,3,2) = ldJ223(i,j,k)
       hes(2,3,3) = ldJ233(i,j,k)

       hes(3,1,1) = ldJ311(i,j,k)
       hes(3,1,2) = ldJ312(i,j,k)
       hes(3,1,3) = ldJ313(i,j,k)
       hes(3,2,1) = ldJ312(i,j,k)
       hes(3,2,2) = ldJ322(i,j,k)
       hes(3,2,3) = ldJ323(i,j,k)
       hes(3,3,1) = ldJ313(i,j,k)
       hes(3,3,2) = ldJ323(i,j,k)
       hes(3,3,3) = ldJ333(i,j,k)
    else
       jac      = 0.0
       jac(1,1) = 1.0
       jac(2,2) = 1.0
       jac(3,3) = 1.0
       hes      = 0.0
    end if
    ! write(*,*) 'J = ', jac

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

    !---- compute beta in local coordinates ----
    beta_l = 0
    do m = 1, 3
       do a = 1, 3
          beta_l(m) = beta_l(m) + jac(m,a) * beta(a)
       end do
    end do
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
    d1_hh11(1) = (   -hxx(i+2,j,k) + 8*hxx(i+1,j,k)                      &
                  - 8*hxx(i-1,j,k) +   hxx(i-2,j,k) ) / dx12
    d1_hh12(1) = (   -hxy(i+2,j,k) + 8*hxy(i+1,j,k)                      &
                  - 8*hxy(i-1,j,k) +   hxy(i-2,j,k) ) / dx12
    d1_hh13(1) = (   -hxz(i+2,j,k) + 8*hxz(i+1,j,k)                      &
                  - 8*hxz(i-1,j,k) +   hxz(i-2,j,k) ) / dx12
    d1_hh22(1) = (   -hyy(i+2,j,k) + 8*hyy(i+1,j,k)                      &
                  - 8*hyy(i-1,j,k) +   hyy(i-2,j,k) ) / dx12
    d1_hh23(1) = (   -hyz(i+2,j,k) + 8*hyz(i+1,j,k)                      &
                  - 8*hyz(i-1,j,k) +   hyz(i-2,j,k) ) / dx12
    d1_hh33(1) = (   -hzz(i+2,j,k) + 8*hzz(i+1,j,k)                      &
                  - 8*hzz(i-1,j,k) +   hzz(i-2,j,k) ) / dx12

    d1_hh11(2) = (   -hxx(i,j+2,k) + 8*hxx(i,j+1,k)                      &
                  - 8*hxx(i,j-1,k) +   hxx(i,j-2,k) ) / dy12
    d1_hh12(2) = (   -hxy(i,j+2,k) + 8*hxy(i,j+1,k)                      &
                  - 8*hxy(i,j-1,k) +   hxy(i,j-2,k) ) / dy12
    d1_hh13(2) = (   -hxz(i,j+2,k) + 8*hxz(i,j+1,k)                      &
                  - 8*hxz(i,j-1,k) +   hxz(i,j-2,k) ) / dy12
    d1_hh22(2) = (   -hyy(i,j+2,k) + 8*hyy(i,j+1,k)                      &
                  - 8*hyy(i,j-1,k) +   hyy(i,j-2,k) ) / dy12
    d1_hh23(2) = (   -hyz(i,j+2,k) + 8*hyz(i,j+1,k)                      &
                  - 8*hyz(i,j-1,k) +   hyz(i,j-2,k) ) / dy12
    d1_hh33(2) = (   -hzz(i,j+2,k) + 8*hzz(i,j+1,k)                      &
                  - 8*hzz(i,j-1,k) +   hzz(i,j-2,k) ) / dy12

    d1_hh11(3) = (   -hxx(i,j,k+2) + 8*hxx(i,j,k+1)                      &
                  - 8*hxx(i,j,k-1) +   hxx(i,j,k-2) ) / dz12
    d1_hh12(3) = (   -hxy(i,j,k+2) + 8*hxy(i,j,k+1)                      &
                  - 8*hxy(i,j,k-1) +   hxy(i,j,k-2) ) / dz12
    d1_hh13(3) = (   -hxz(i,j,k+2) + 8*hxz(i,j,k+1)                      &
                  - 8*hxz(i,j,k-1) +   hxz(i,j,k-2) ) / dz12
    d1_hh22(3) = (   -hyy(i,j,k+2) + 8*hyy(i,j,k+1)                      &
                  - 8*hyy(i,j,k-1) +   hyy(i,j,k-2) ) / dz12
    d1_hh23(3) = (   -hyz(i,j,k+2) + 8*hyz(i,j,k+1)                      &
                  - 8*hyz(i,j,k-1) +   hyz(i,j,k-2) ) / dz12
    d1_hh33(3) = (   -hzz(i,j,k+2) + 8*hzz(i,j,k+1)                      &
                  - 8*hzz(i,j,k-1) +   hzz(i,j,k-2) ) / dz12


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

      di = int( sign( one, beta_l(1) ) )
      dj = int( sign( one, beta_l(2) ) )
      dk = int( sign( one, beta_l(3) ) )

      ! ad1_lphi1
      d1_f(1) = di * ( -3*phi1(i-di,j,k) - 10*phi1(i,j,k) + 18*phi1(i+di,j,k)   &
                      - 6*phi1(i+2*di,j,k)  + phi1(i+3*di,j,k)) / dx12
      d1_f(2) = dj * ( -3*phi1(i,j-dj,k) - 10*phi1(i,j,k) + 18*phi1(i,j+dj,k)   &
                      - 6*phi1(i,j+2*dj,k)  + phi1(i,j+3*dj,k)) / dy12
      d1_f(3) = dk * ( -3*phi1(i,j,k-dk) - 10*phi1(i,j,k) + 18*phi1(i,j,k+dk)   &
                      - 6*phi1(i,j,k+2*dk)  + phi1(i,j,k+3*dk)) / dz12
      ad1_lphi1 = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

      ! ad1_lKphi1
      d1_f(1) = di * ( -3*Kphi1(i-di,j,k) - 10*Kphi1(i,j,k) + 18*Kphi1(i+di,j,k)   &
                      - 6*Kphi1(i+2*di,j,k)  + Kphi1(i+3*di,j,k)) / dx12
      d1_f(2) = dj * ( -3*Kphi1(i,j-dj,k) - 10*Kphi1(i,j,k) + 18*Kphi1(i,j+dj,k)   &
                      - 6*Kphi1(i,j+2*dj,k)  + Kphi1(i,j+3*dj,k)) / dy12
      d1_f(3) = dk * ( -3*Kphi1(i,j,k-dk) - 10*Kphi1(i,j,k) + 18*Kphi1(i,j,k+dk)   &
                      - 6*Kphi1(i,j,k+2*dk)  + Kphi1(i,j,k+3*dk)) / dz12
      ad1_lKphi1 = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

      ! ad1_lphi2
      d1_f(1) = di * ( -3*phi2(i-di,j,k) - 10*phi2(i,j,k) + 18*phi2(i+di,j,k)   &
                      - 6*phi2(i+2*di,j,k)  + phi2(i+3*di,j,k)) / dx12
      d1_f(2) = dj * ( -3*phi2(i,j-dj,k) - 10*phi2(i,j,k) + 18*phi2(i,j+dj,k)   &
                      - 6*phi2(i,j+2*dj,k)  + phi2(i,j+3*dj,k)) / dy12
      d1_f(3) = dk * ( -3*phi2(i,j,k-dk) - 10*phi2(i,j,k) + 18*phi2(i,j,k+dk)   &
                      - 6*phi2(i,j,k+2*dk)  + phi2(i,j,k+3*dk)) / dz12
      ad1_lphi2 = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

      ! ad1_lKphi2
      d1_f(1) = di * ( -3*Kphi2(i-di,j,k) - 10*Kphi2(i,j,k) + 18*Kphi2(i+di,j,k)   &
                      - 6*Kphi2(i+2*di,j,k)  + Kphi2(i+3*di,j,k)) / dx12
      d1_f(2) = dj * ( -3*Kphi2(i,j-dj,k) - 10*Kphi2(i,j,k) + 18*Kphi2(i,j+dj,k)   &
                      - 6*Kphi2(i,j+2*dj,k)  + Kphi2(i,j+3*dj,k)) / dy12
      d1_f(3) = dk * ( -3*Kphi2(i,j,k-dk) - 10*Kphi2(i,j,k) + 18*Kphi2(i,j,k+dk)   &
                      - 6*Kphi2(i,j,k+2*dk)  + Kphi2(i,j,k+3*dk)) / dz12
      ad1_lKphi2 = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

    else

      ! ad1_lphi1
      ad1_lphi1 = beta_l(1)*d1_lphi1(1) + beta_l(2)*d1_lphi1(2) + beta_l(3)*d1_lphi1(3)

      ! ad1_lKphi1
      ad1_lKphi1 = beta_l(1)*d1_lKphi1(1) + beta_l(2)*d1_lKphi1(2) + beta_l(3)*d1_lKphi1(3)

      ! ad1_lphi2
      ad1_lphi2 = beta_l(1)*d1_lphi2(1) + beta_l(2)*d1_lphi2(2) + beta_l(3)*d1_lphi2(3)

      ! ad1_lKphi2
      ad1_lKphi2 = beta_l(1)*d1_lKphi2(1) + beta_l(2)*d1_lKphi2(2) + beta_l(3)*d1_lKphi2(3)

    end if

    !-------------------------------------------
    if (use_jacobian) then
       call Scalar_apply_jacobian(d1_alph, jac)
       call Scalar_apply_jacobian(d1_lKphi1, jac)
       call Scalar_apply_jacobian(d1_lKphi2, jac)
       call Scalar_apply_jacobian(d1_hh11, jac)
       call Scalar_apply_jacobian(d1_hh12, jac)
       call Scalar_apply_jacobian(d1_hh13, jac)
       call Scalar_apply_jacobian(d1_hh22, jac)
       call Scalar_apply_jacobian(d1_hh23, jac)
       call Scalar_apply_jacobian(d1_hh33, jac)

       call Scalar_apply_jacobian2(d1_lphi1, d2_lphi1, jac, hes)
       call Scalar_apply_jacobian2(d1_lphi2, d2_lphi2, jac, hes)
    end if

    d1_hh(1,1,:) = d1_hh11(:)
    d1_hh(1,2,:) = d1_hh12(:)
    d1_hh(1,3,:) = d1_hh13(:)
    d1_hh(2,2,:) = d1_hh22(:)
    d1_hh(2,3,:) = d1_hh23(:)
    d1_hh(3,3,:) = d1_hh33(:)
    d1_hh(2,1,:) = d1_hh(1,2,:)
    d1_hh(3,1,:) = d1_hh(1,3,:)
    d1_hh(3,2,:) = d1_hh(2,3,:)
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
!
!=============================================================================
!
subroutine Scalar_ord4_calc_rhs_bdry_sph( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT  i, j, k
  CCTK_REAL odr2
  CCTK_REAL alph, lphi1, lphi2, lKphi1, lKphi2
  CCTK_REAL dr_lphi1, dr_lphi2, dr_lKphi1, dr_lKphi2
  CCTK_REAL rr

  CCTK_INT  reflevel, map

  odr2 = 1 / (2*CCTK_DELTA_SPACE(3))

  reflevel = GetRefinementLevel(cctkGH)
  map      = MultiPatch_GetMap(cctkGH)

  ! apply only on the coarsest level and in the spherical shell. points marked
  ! with cctk_bbox(6) == 0 are inter-processor boundaries, so we also do not
  ! want those
  if (reflevel /= 0 .or. map == 0 .or. cctk_bbox(6) == 0) return

  do k = cctk_lsh(3)-cctk_nghostzones(3)+1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

           rr       = r(i,j,k)

           alph     = alp(i,j,k)
           lphi1    = phi1(i,j,k)
           lphi2    = phi2(i,j,k)
           lKphi1   = Kphi1(i,j,k)
           lKphi2   = Kphi2(i,j,k)

           dr_lphi1  = (phi1(i,j,k-2) - 4*phi1(i,j,k-1) + 3*phi1(i,j,k))*odr2
           dr_lphi2  = (phi2(i,j,k-2) - 4*phi2(i,j,k-1) + 3*phi2(i,j,k))*odr2
           dr_lKphi1 = (Kphi1(i,j,k-2) - 4*Kphi1(i,j,k-1) + 3*Kphi1(i,j,k))*odr2
           dr_lKphi2 = (Kphi2(i,j,k-2) - 4*Kphi2(i,j,k-1) + 3*Kphi2(i,j,k))*odr2

           ! FIXME: the boundary conditions here should really be the ones
           ! below, but more tests are needed
           rhs_phi1(i,j,k)  = -dr_lphi1 - lphi1 / rr
           rhs_phi2(i,j,k)  = -dr_lphi2 - lphi2 / rr
           ! rhs_phi1(i,j,k)  = -2.0 * alph * lKphi1
           ! rhs_phi2(i,j,k)  = -2.0 * alph * lKphi2

           rhs_Kphi1(i,j,k) = -dr_lKphi1 - lKphi1 / rr + 0.25 * mu**2 * alph * lphi1
           rhs_Kphi2(i,j,k) = -dr_lKphi2 - lKphi2 / rr + 0.25 * mu**2 * alph * lphi2

        end do
     end do
  end do

end subroutine Scalar_ord4_calc_rhs_bdry_sph
