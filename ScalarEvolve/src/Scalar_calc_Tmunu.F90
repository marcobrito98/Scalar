#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

subroutine Scalar_calc_Tmunu( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  ! Fundamental variables
  CCTK_REAL                alph, beta(3), betad(3)
  CCTK_REAL                gg(4,4), gu(4,4), deth
  CCTK_REAL                lphi1, lphi2, lKphi1, lKphi2
  CCTK_COMPLEX             lphi, absphi2
  CCTK_REAL                Tab(4,4)

  ! Fluxes variables
  CCTK_REAL                rhoSF, jrSF, jiSF(3), SrrSF, SijSF(3,3)
  CCTK_REAL                xx(3), rr

  ! First derivatives
  CCTK_REAL                d1_lphi1(4), d1_lphi2(4)
  CCTK_COMPLEX             d1_lphi(4)

  ! Misc variables
  CCTK_REAL                dx12, dy12, dz12
  CCTK_REAL                aux

  CCTK_INT                 i, j, k
  CCTK_INT                 a, b, c

  ! jacobian
  integer                  istat
  logical                  use_jacobian
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ11, lJ12, lJ13
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ21, lJ22, lJ23
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ31, lJ32, lJ33
  CCTK_POINTER             lJ11_ptr, lJ12_ptr, lJ13_ptr
  CCTK_POINTER             lJ21_ptr, lJ22_ptr, lJ23_ptr
  CCTK_POINTER             lJ31_ptr, lJ32_ptr, lJ33_ptr
  CCTK_REAL                jac(3,3)

  pointer (lJ11_ptr, lJ11), (lJ12_ptr, lJ12), (lJ13_ptr, lJ13)
  pointer (lJ21_ptr, lJ21), (lJ22_ptr, lJ22), (lJ23_ptr, lJ23)
  pointer (lJ31_ptr, lJ31), (lJ32_ptr, lJ32), (lJ33_ptr, lJ33)

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
  end if

  dx12 = 12*CCTK_DELTA_SPACE(1)
  dy12 = 12*CCTK_DELTA_SPACE(2)
  dz12 = 12*CCTK_DELTA_SPACE(3)


  !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(k,j,i, a,b, aux,&
  !$OMP                                 alph,beta,betad,&
  !$OMP                                 gg,gu,deth,&
  !$OMP                                 lphi1,lphi2,lKphi1,lKphi2,&
  !$OMP                                 lphi,absphi2,d1_lphi,&
  !$OMP                                 rhoSF, jrSF, jiSF, SrrSF, SijSF,&
  !$OMP                                 Tab, d1_lphi1, d1_lphi2, jac)
  do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
     do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
        do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)

           !------------ Get local variables ----------
           lphi1     = phi1(i,j,k)
           lphi2     = phi2(i,j,k)
           lKphi1    = Kphi1(i,j,k)
           lKphi2    = Kphi2(i,j,k)

           alph      = alp(i,j,k)

           beta(1)   = betax(i,j,k)
           beta(2)   = betay(i,j,k)
           beta(3)   = betaz(i,j,k)

           gg(1,1)   = gxx(i,j,k)
           gg(1,2)   = gxy(i,j,k)
           gg(1,3)   = gxz(i,j,k)
           gg(2,2)   = gyy(i,j,k)
           gg(2,3)   = gyz(i,j,k)
           gg(3,3)   = gzz(i,j,k)
           gg(2,1)   = gg(1,2)
           gg(3,1)   = gg(1,3)
           gg(3,2)   = gg(2,3)

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
           else
              jac      = 0.0d0
              jac(1,1) = 1.0d0
              jac(2,2) = 1.0d0
              jac(3,3) = 1.0d0
           end if

           ! now we compute beta_i (betad)
           betad = 0.0d0
           do a = 1, 3
              do b = 1, 3
                 betad(a) = betad(a) + gg(a,b) * beta(b)
              end do
           end do

           ! and finish with the rest of the 4-metric gg.
           ! we will use the 4th slot for time throughout
           gg(1,4)   = betad(1)
           gg(2,4)   = betad(2)
           gg(3,4)   = betad(3)
           gg(4,1)   = gg(1,4)
           gg(4,2)   = gg(2,4)
           gg(4,3)   = gg(3,4)

           gg(4,4)   = -alph * alph
           do a = 1, 3
              gg(4,4) = gg(4,4) + beta(a) * betad(a)
           end do

           !------------ Invert metric ----------------

           ! determinant of the 3-metric
           deth    =     gg(1,1) * gg(2,2) * gg(3,3)                              &
                   + 2 * gg(1,2) * gg(1,3) * gg(2,3)                              &
                   -     gg(1,1) * gg(2,3) ** 2                                   &
                   -     gg(2,2) * gg(1,3) ** 2                                   &
                   -     gg(3,3) * gg(1,2) ** 2

           gu(1,1) = (gg(2,2) * gg(3,3) - gg(2,3) ** 2     ) / deth               &
                   - beta(1) * beta(1) / (alph * alph)
           gu(2,2) = (gg(1,1) * gg(3,3) - gg(1,3) ** 2     ) / deth               &
                   - beta(2) * beta(2) / (alph * alph)
           gu(3,3) = (gg(1,1) * gg(2,2) - gg(1,2) ** 2     ) / deth               &
                   - beta(3) * beta(3) / (alph * alph)
           gu(1,2) = (gg(1,3) * gg(2,3) - gg(1,2) * gg(3,3)) / deth               &
                   - beta(1) * beta(2) / (alph * alph)
           gu(1,3) = (gg(1,2) * gg(2,3) - gg(1,3) * gg(2,2)) / deth               &
                   - beta(1) * beta(3) / (alph * alph)
           gu(2,3) = (gg(1,3) * gg(1,2) - gg(2,3) * gg(1,1)) / deth               &
                   - beta(2) * beta(3) / (alph * alph)

           gu(1,4) = beta(1) / (alph * alph)
           gu(2,4) = beta(2) / (alph * alph)
           gu(3,4) = beta(3) / (alph * alph)
           gu(4,4) = -1.0d0 / (alph * alph)

           gu(2,1) = gu(1,2)
           gu(3,1) = gu(1,3)
           gu(3,2) = gu(2,3)

           gu(4,1) = gu(1,4)
           gu(4,2) = gu(2,4)
           gu(4,3) = gu(3,4)


           !------------- Centered 1st derivatives -----------

           ! d1_lphi1(3)
           d1_lphi1(1)  = (   -phi1(i+2,j,k) + 8*phi1(i+1,j,k)                        &
                           - 8*phi1(i-1,j,k) +   phi1(i-2,j,k) ) / dx12

           d1_lphi1(2)  = (   -phi1(i,j+2,k) + 8*phi1(i,j+1,k)                        &
                           - 8*phi1(i,j-1,k) +   phi1(i,j-2,k) ) / dy12

           d1_lphi1(3)  = (   -phi1(i,j,k+2) + 8*phi1(i,j,k+1)                        &
                           - 8*phi1(i,j,k-1) +   phi1(i,j,k-2) ) / dz12

           ! d1_lphi2(3)
           d1_lphi2(1)  = (   -phi2(i+2,j,k) + 8*phi2(i+1,j,k)                        &
                           - 8*phi2(i-1,j,k) +   phi2(i-2,j,k) ) / dx12

           d1_lphi2(2)  = (   -phi2(i,j+2,k) + 8*phi2(i,j+1,k)                        &
                           - 8*phi2(i,j-1,k) +   phi2(i,j-2,k) ) / dy12

           d1_lphi2(3)  = (   -phi2(i,j,k+2) + 8*phi2(i,j,k+1)                        &
                           - 8*phi2(i,j,k-1) +   phi2(i,j,k-2) ) / dz12

           !-------------------------------------------
           if (use_jacobian) then
              call Scalar_apply_jacobian(d1_lphi1, jac)
              call Scalar_apply_jacobian(d1_lphi2, jac)
           end if
           !-------------------------------------------

           ! time derivatives
           d1_lphi1(4)  = -2 * alph * lKphi1
           d1_lphi2(4)  = -2 * alph * lKphi2
           do a = 1, 3
              d1_lphi1(4) = d1_lphi1(4) + beta(a) * d1_lphi1(a)
              d1_lphi2(4) = d1_lphi2(4) + beta(a) * d1_lphi2(a)
           end do

           !-------------------------------------------

           ! build complex scalar field
           lphi     = dcmplx(lphi1, lphi2)
           d1_lphi  = dcmplx(d1_lphi1, d1_lphi2)
           absphi2  = real( lphi * conjg(lphi) )

           aux = mu * mu * absphi2                                            &
                * (1 - 2 * V_lambda * absphi2) * (1 - 2 * V_lambda * absphi2)

           do a = 1, 4
              do b = 1, 4
                 aux = aux + gu(a,b) * real( d1_lphi(a) * conjg(d1_lphi(b)) )
              end do
           end do

           ! compute the stress-energy tensor
           do a = 1, 4
              do b = 1, 4
                 Tab(a,b) = real(        d1_lphi(a)  * conjg(d1_lphi(b))     &
                                 + conjg(d1_lphi(a)) *       d1_lphi(b) )    &
                          - gg(a,b) * aux
              end do
           end do

           ! and finally store it in the Tmunu variables
           eTtt(i,j,k) = eTtt(i,j,k) + Tab(4,4)
           eTtx(i,j,k) = eTtx(i,j,k) + Tab(4,1)
           eTty(i,j,k) = eTty(i,j,k) + Tab(4,2)
           eTtz(i,j,k) = eTtz(i,j,k) + Tab(4,3)
           eTxx(i,j,k) = eTxx(i,j,k) + Tab(1,1)
           eTxy(i,j,k) = eTxy(i,j,k) + Tab(1,2)
           eTxz(i,j,k) = eTxz(i,j,k) + Tab(1,3)
           eTyy(i,j,k) = eTyy(i,j,k) + Tab(2,2)
           eTyz(i,j,k) = eTyz(i,j,k) + Tab(2,3)
           eTzz(i,j,k) = eTzz(i,j,k) + Tab(3,3)

           if ( ( compute_fluxes == 1 ) .and. ( use_jacobian .eqv. .false. ) ) then
                !--- local coordinates to define local radial coordinate --
                xx(1) = x(i,j,k)
                xx(2) = y(i,j,k)
                xx(3) = z(i,j,k)

                rr = sqrt( xx(1)**2 + xx(2)**2 + xx(3)**2 )
                if( rr < eps_r ) rr = eps_r
                !----------------------------------------------------------

                !--- Eulerian energy density ------------------------------
                rhoSF = 0
                rhoSF = Tab(4,4) / ( alph * alph )

                !--- Eulerian energy-momentum flux ------------------------
                ! jrSF
                ! calculate jiSF(3) in Cartesian coordinates

                jiSF = 0
                do a = 1, 3
                    do b = 1, 3
                        do c = 1, 3
                            jiSF(a) = jiSF(a) + gg(a,b) * gu(b,c) * Tab(c,4)
                        end do
                    end do
                end do
                !jiSF =  - jiSF
                jiSF =  - jiSF / alph

                !--- Trasformation to radial spherical coordinate ---------
                jrSF = 0
                do a = 1, 3
                    jrSF = jrSF + jiSF(a) * xx(a) / rr
                end do
                !----------------------------------------------------------
                !--- Spatial stress tensor --------------------------------
                ! SrrSF
                ! calculate SijSF(3,3) in Cartesian coordinates

                !SijSF = 0
                !do a = 1, 3
                !  do b = 1, 3
                !    SijSF(a,b) = SijSF(a,b) + Tab(a,b)
                !  end do
                !end do

                !--- Trasformation to radial spherical coordinate ---------
                SrrSF = 0
                do a = 1, 3
                  do b = 1, 3
                    SrrSF = SrrSF + Tab(a,b) * xx(a) / rr * xx(b) / rr
                  end do
                end do
                !----------------------------------------------------------

                !------------ Write to grid functions ---------------------
                rhoSF_gf(i,j,k)    = rhoSF
                jrSF_gf(i,j,k)     = jrSF
                SrrSF_gf(i,j,k)    = SrrSF
                !----------------------------------------------------------

           end if

        end do
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine Scalar_calc_Tmunu
