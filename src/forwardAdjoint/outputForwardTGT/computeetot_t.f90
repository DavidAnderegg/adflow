   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.7 (r4786) - 21 Feb 2013 15:53
   !
   !  Differentiation of computeetot in forward (tangent) mode (with options debugTangent i4 dr8 r8):
   !   variations   of useful results: *gamma *w
   !   with respect to varying inputs: *p *gamma *w rgas
   !   Plus diff mem management of: p:in gamma:in w:in
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          computeEtot.F90                                 *
   !      * Author:        Edwin van der Weide, Steve Repsher              *
   !      * Starting date: 08-13-2003                                      *
   !      * Last modified: 10-14-2005                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE COMPUTEETOT_T(istart, iend, jstart, jend, kstart, kend, &
   &  correctfork)
   USE FLOWVARREFSTATE
   USE BLOCKPOINTERS_D
   USE INPUTPHYSICS
   USE DIFFSIZES
   !  Hint: ISIZE4OFDrfw should be the size of dimension 4 of array *w
   !  Hint: ISIZE3OFDrfw should be the size of dimension 3 of array *w
   !  Hint: ISIZE2OFDrfw should be the size of dimension 2 of array *w
   !  Hint: ISIZE1OFDrfw should be the size of dimension 1 of array *w
   !  Hint: ISIZE3OFDrfp should be the size of dimension 3 of array *p
   !  Hint: ISIZE2OFDrfp should be the size of dimension 2 of array *p
   !  Hint: ISIZE1OFDrfp should be the size of dimension 1 of array *p
   !  Hint: ISIZE3OFDrfgamma should be the size of dimension 3 of array *gamma
   !  Hint: ISIZE2OFDrfgamma should be the size of dimension 2 of array *gamma
   !  Hint: ISIZE1OFDrfgamma should be the size of dimension 1 of array *gamma
   IMPLICIT NONE
   !
   !      ******************************************************************
   !      *                                                                *
   !      * ComputeEtot computes the total energy from the given density,  *
   !      * velocity and presssure. For a calorically and thermally        *
   !      * perfect gas the well-known expression is used; for only a      *
   !      * thermally perfect gas, cp is a function of temperature, curve  *
   !      * fits are used and a more complex expression is obtained.       *
   !      * It is assumed that the pointers in blockPointers already       *
   !      * point to the correct block.                                    *
   !      *                                                                *
   !      ******************************************************************
   !
   !
   !      Subroutine arguments.
   !
   INTEGER(kind=inttype), INTENT(IN) :: istart, iend, jstart, jend
   INTEGER(kind=inttype), INTENT(IN) :: kstart, kend
   LOGICAL, INTENT(IN) :: correctfork
   !
   !      Local variables.
   !
   INTEGER(kind=inttype) :: i, j, k
   REAL(kind=realtype) :: ovgm1, factk, scale
   REAL(kind=realtype) :: scaled
   EXTERNAL DEBUG_TGT_HERE
   LOGICAL :: DEBUG_TGT_HERE
   IF (.TRUE. .AND. DEBUG_TGT_HERE('entry', .FALSE.)) THEN
   CALL DEBUG_TGT_REAL8ARRAY('p', p, pd, ISIZE1OFDrfp*ISIZE2OFDrfp*&
   &                        ISIZE3OFDrfp)
   CALL DEBUG_TGT_REAL8ARRAY('w', w, wd, ISIZE1OFDrfw*ISIZE2OFDrfw*&
   &                        ISIZE3OFDrfw*ISIZE4OFDrfw)
   CALL DEBUG_TGT_REAL8('rgas', rgas, rgasd)
   CALL DEBUG_TGT_DISPLAY('entry')
   END IF
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Determine the cp model used in the computation.
   SELECT CASE  (cpmodel) 
   CASE (cpconstant) 
   ! Constant cp and thus constant gamma.
   ! Abbreviate 1/(gamma -1) a bit easier.
   ovgm1 = one/(gammaconstant-one)
   ! Loop over the given range of the block and compute the first
   ! step of the energy.
   DO k=kstart,kend
   DO j=jstart,jend
   DO i=istart,iend
   wd(i, j, k, irhoe) = ovgm1*pd(i, j, k) + half*(wd(i, j, k, &
   &            irho)*(w(i, j, k, ivx)**2+w(i, j, k, ivy)**2+w(i, j, k, ivz)&
   &            **2)+w(i, j, k, irho)*(2*w(i, j, k, ivx)*wd(i, j, k, ivx)+2*&
   &            w(i, j, k, ivy)*wd(i, j, k, ivy)+2*w(i, j, k, ivz)*wd(i, j, &
   &            k, ivz)))
   w(i, j, k, irhoe) = ovgm1*p(i, j, k) + half*w(i, j, k, irho)*(&
   &            w(i, j, k, ivx)**2+w(i, j, k, ivy)**2+w(i, j, k, ivz)**2)
   END DO
   END DO
   END DO
   ! Second step. Correct the energy in case a turbulent kinetic
   ! energy is present.
   IF (correctfork) THEN
   IF (.TRUE. .AND. DEBUG_TGT_HERE('middle', .FALSE.)) THEN
   CALL DEBUG_TGT_REAL8ARRAY('gamma', gamma, gammad, &
   &                            ISIZE1OFDrfgamma*ISIZE2OFDrfgamma*&
   &                            ISIZE3OFDrfgamma)
   CALL DEBUG_TGT_REAL8ARRAY('w', w, wd, ISIZE1OFDrfw*ISIZE2OFDrfw*&
   &                            ISIZE3OFDrfw*ISIZE4OFDrfw)
   CALL DEBUG_TGT_DISPLAY('middle')
   END IF
   factk = ovgm1*(five*third-gammaconstant)
   DO k=kstart,kend
   DO j=jstart,jend
   DO i=istart,iend
   wd(i, j, k, irhoe) = wd(i, j, k, irhoe) - factk*(wd(i, j, k&
   &              , irho)*w(i, j, k, itu1)+w(i, j, k, irho)*wd(i, j, k, itu1&
   &              ))
   w(i, j, k, irhoe) = w(i, j, k, irhoe) - factk*w(i, j, k, &
   &              irho)*w(i, j, k, itu1)
   END DO
   END DO
   END DO
   END IF
   CASE (cptempcurvefits) 
   !        ================================================================
   ! Cp as function of the temperature is given via curve fits.
   ! Store a scale factor to compute the nonDimensional
   ! internal energy.
   scaled = rgasd/tref
   scale = rgas/tref
   ! Loop over the given range of the block.
   DO k=kstart,kend
   DO j=jstart,jend
   DO i=istart,iend
   CALL DEBUG_TGT_CALL('COMPUTEETOTCELLCPFIT', .TRUE., .FALSE.)
   CALL COMPUTEETOTCELLCPFIT_T(i, j, k, scale, scaled, &
   &                                correctfork)
   CALL DEBUG_TGT_EXIT()
   END DO
   END DO
   END DO
   END SELECT
   IF (.TRUE. .AND. DEBUG_TGT_HERE('exit', .FALSE.)) THEN
   CALL DEBUG_TGT_REAL8ARRAY('gamma', gamma, gammad, ISIZE1OFDrfgamma*&
   &                        ISIZE2OFDrfgamma*ISIZE3OFDrfgamma)
   CALL DEBUG_TGT_REAL8ARRAY('w', w, wd, ISIZE1OFDrfw*ISIZE2OFDrfw*&
   &                        ISIZE3OFDrfw*ISIZE4OFDrfw)
   CALL DEBUG_TGT_DISPLAY('exit')
   END IF
   !
   40 FORMAT(1x,i4,i4,i4,e20.6)
   END SUBROUTINE COMPUTEETOT_T
