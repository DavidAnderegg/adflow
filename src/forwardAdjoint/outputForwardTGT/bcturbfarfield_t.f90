   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.7 (r4786) - 21 Feb 2013 15:53
   !
   !  Differentiation of bcturbfarfield in forward (tangent) mode (with options debugTangent i4 dr8 r8):
   !   variations   of useful results: *bvtj1 *bvtj2 *bmtk1 *bmtk2
   !                *bvtk1 *bvtk2 *bmti1 *bmti2 *bvti1 *bvti2 *bmtj1
   !                *bmtj2
   !   with respect to varying inputs: *bvtj1 *bvtj2 *bmtk1 *bmtk2
   !                *bvtk1 *bvtk2 *bmti1 *bmti2 *bvti1 *bvti2 *bmtj1
   !                *bmtj2 winf
   !   Plus diff mem management of: bvtj1:in bvtj2:in bmtk1:in bmtk2:in
   !                bvtk1:in bvtk2:in bmti1:in bmti2:in bvti1:in bvti2:in
   !                bmtj1:in bmtj2:in bcdata:in *bcdata.norm:in
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          bcTurbFarfield.f90                              *
   !      * Author:        Georgi Kalitzin, Edwin van der Weide            *
   !      * Starting date: 06-15-2003                                      *
   !      * Last modified: 06-12-2005                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE BCTURBFARFIELD_T(nn)
   USE CONSTANTS
   USE FLOWVARREFSTATE
   USE BLOCKPOINTERS_D
   USE BCTYPES
   USE DIFFSIZES
   !  Hint: ISIZE4OFDrfbmtj2 should be the size of dimension 4 of array *bmtj2
   !  Hint: ISIZE3OFDrfbmtj2 should be the size of dimension 3 of array *bmtj2
   !  Hint: ISIZE2OFDrfbmtj2 should be the size of dimension 2 of array *bmtj2
   !  Hint: ISIZE1OFDrfbmtj2 should be the size of dimension 1 of array *bmtj2
   !  Hint: ISIZE4OFDrfbmtj1 should be the size of dimension 4 of array *bmtj1
   !  Hint: ISIZE3OFDrfbmtj1 should be the size of dimension 3 of array *bmtj1
   !  Hint: ISIZE2OFDrfbmtj1 should be the size of dimension 2 of array *bmtj1
   !  Hint: ISIZE1OFDrfbmtj1 should be the size of dimension 1 of array *bmtj1
   !  Hint: ISIZE3OFDrfbvti2 should be the size of dimension 3 of array *bvti2
   !  Hint: ISIZE2OFDrfbvti2 should be the size of dimension 2 of array *bvti2
   !  Hint: ISIZE1OFDrfbvti2 should be the size of dimension 1 of array *bvti2
   !  Hint: ISIZE3OFDrfbvti1 should be the size of dimension 3 of array *bvti1
   !  Hint: ISIZE2OFDrfbvti1 should be the size of dimension 2 of array *bvti1
   !  Hint: ISIZE1OFDrfbvti1 should be the size of dimension 1 of array *bvti1
   !  Hint: ISIZE4OFDrfbmti2 should be the size of dimension 4 of array *bmti2
   !  Hint: ISIZE3OFDrfbmti2 should be the size of dimension 3 of array *bmti2
   !  Hint: ISIZE2OFDrfbmti2 should be the size of dimension 2 of array *bmti2
   !  Hint: ISIZE1OFDrfbmti2 should be the size of dimension 1 of array *bmti2
   !  Hint: ISIZE4OFDrfbmti1 should be the size of dimension 4 of array *bmti1
   !  Hint: ISIZE3OFDrfbmti1 should be the size of dimension 3 of array *bmti1
   !  Hint: ISIZE2OFDrfbmti1 should be the size of dimension 2 of array *bmti1
   !  Hint: ISIZE1OFDrfbmti1 should be the size of dimension 1 of array *bmti1
   !  Hint: ISIZE3OFDrfbvtk2 should be the size of dimension 3 of array *bvtk2
   !  Hint: ISIZE2OFDrfbvtk2 should be the size of dimension 2 of array *bvtk2
   !  Hint: ISIZE1OFDrfbvtk2 should be the size of dimension 1 of array *bvtk2
   !  Hint: ISIZE3OFDrfbvtk1 should be the size of dimension 3 of array *bvtk1
   !  Hint: ISIZE2OFDrfbvtk1 should be the size of dimension 2 of array *bvtk1
   !  Hint: ISIZE1OFDrfbvtk1 should be the size of dimension 1 of array *bvtk1
   !  Hint: ISIZE4OFDrfbmtk2 should be the size of dimension 4 of array *bmtk2
   !  Hint: ISIZE3OFDrfbmtk2 should be the size of dimension 3 of array *bmtk2
   !  Hint: ISIZE2OFDrfbmtk2 should be the size of dimension 2 of array *bmtk2
   !  Hint: ISIZE1OFDrfbmtk2 should be the size of dimension 1 of array *bmtk2
   !  Hint: ISIZE4OFDrfbmtk1 should be the size of dimension 4 of array *bmtk1
   !  Hint: ISIZE3OFDrfbmtk1 should be the size of dimension 3 of array *bmtk1
   !  Hint: ISIZE2OFDrfbmtk1 should be the size of dimension 2 of array *bmtk1
   !  Hint: ISIZE1OFDrfbmtk1 should be the size of dimension 1 of array *bmtk1
   !  Hint: ISIZE3OFDrfbvtj2 should be the size of dimension 3 of array *bvtj2
   !  Hint: ISIZE2OFDrfbvtj2 should be the size of dimension 2 of array *bvtj2
   !  Hint: ISIZE1OFDrfbvtj2 should be the size of dimension 1 of array *bvtj2
   !  Hint: ISIZE3OFDrfbvtj1 should be the size of dimension 3 of array *bvtj1
   !  Hint: ISIZE2OFDrfbvtj1 should be the size of dimension 2 of array *bvtj1
   !  Hint: ISIZE1OFDrfbvtj1 should be the size of dimension 1 of array *bvtj1
   IMPLICIT NONE
   !
   !      ******************************************************************
   !      *                                                                *
   !      * bcTurbFarfield applies the implicit treatment of the           *
   !      * farfield boundary condition to subface nn. As the farfield     *
   !      * boundary condition is independent of the turbulence model,     *
   !      * this routine is valid for all models. It is assumed that the   *
   !      * pointers in blockPointers are already set to the correct       *
   !      * block on the correct grid level.                               *
   !      *                                                                *
   !      ******************************************************************
   !
   !
   !      Subroutine arguments.
   !
   INTEGER(kind=inttype), INTENT(IN) :: nn
   !
   !      Local variables.
   !
   INTEGER(kind=inttype) :: i, j, l
   REAL(kind=realtype) :: nnx, nny, nnz, dot
   REAL(kind=realtype), DIMENSION(:, :, :, :), POINTER :: bmt
   REAL(kind=realtype), DIMENSION(:, :, :, :), POINTER :: bmtd
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: bvt
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: bvtd
   EXTERNAL DEBUG_TGT_HERE
   LOGICAL :: DEBUG_TGT_HERE
   IF (.TRUE. .AND. DEBUG_TGT_HERE('entry', .FALSE.)) THEN
   CALL DEBUG_TGT_REAL8ARRAY('winf', winf, winfd, 10)
   CALL DEBUG_TGT_DISPLAY('entry')
   END IF
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Set the pointers for bmt and bvt, depending on the block face
   ! on which the subface is located.
   SELECT CASE  (bcfaceid(nn)) 
   CASE (imin) 
   bmtd => bmti1d
   bmt => bmti1
   bvtd => bvti1d
   bvt => bvti1
   CASE (imax) 
   bmtd => bmti2d
   bmt => bmti2
   bvtd => bvti2d
   bvt => bvti2
   CASE (jmin) 
   bmtd => bmtj1d
   bmt => bmtj1
   bvtd => bvtj1d
   bvt => bvtj1
   CASE (jmax) 
   bmtd => bmtj2d
   bmt => bmtj2
   bvtd => bvtj2d
   bvt => bvtj2
   CASE (kmin) 
   bmtd => bmtk1d
   bmt => bmtk1
   bvtd => bvtk1d
   bvt => bvtk1
   CASE (kmax) 
   IF (.TRUE. .AND. DEBUG_TGT_HERE('middle', .FALSE.)) THEN
   CALL DEBUG_TGT_REAL8ARRAY('bvtj1', bvtj1, bvtj1d, ISIZE1OFDrfbvtj1&
   &                          *ISIZE2OFDrfbvtj1*ISIZE3OFDrfbvtj1)
   CALL DEBUG_TGT_REAL8ARRAY('bvtj2', bvtj2, bvtj2d, ISIZE1OFDrfbvtj2&
   &                          *ISIZE2OFDrfbvtj2*ISIZE3OFDrfbvtj2)
   CALL DEBUG_TGT_REAL8ARRAY('bmtk1', bmtk1, bmtk1d, ISIZE1OFDrfbmtk1&
   &                          *ISIZE2OFDrfbmtk1*ISIZE3OFDrfbmtk1*&
   &                          ISIZE4OFDrfbmtk1)
   CALL DEBUG_TGT_REAL8ARRAY('bmtk2', bmtk2, bmtk2d, ISIZE1OFDrfbmtk2&
   &                          *ISIZE2OFDrfbmtk2*ISIZE3OFDrfbmtk2*&
   &                          ISIZE4OFDrfbmtk2)
   CALL DEBUG_TGT_REAL8ARRAY('bvtk1', bvtk1, bvtk1d, ISIZE1OFDrfbvtk1&
   &                          *ISIZE2OFDrfbvtk1*ISIZE3OFDrfbvtk1)
   CALL DEBUG_TGT_REAL8ARRAY('bvtk2', bvtk2, bvtk2d, ISIZE1OFDrfbvtk2&
   &                          *ISIZE2OFDrfbvtk2*ISIZE3OFDrfbvtk2)
   CALL DEBUG_TGT_REAL8ARRAY('bmti1', bmti1, bmti1d, ISIZE1OFDrfbmti1&
   &                          *ISIZE2OFDrfbmti1*ISIZE3OFDrfbmti1*&
   &                          ISIZE4OFDrfbmti1)
   CALL DEBUG_TGT_REAL8ARRAY('bmti2', bmti2, bmti2d, ISIZE1OFDrfbmti2&
   &                          *ISIZE2OFDrfbmti2*ISIZE3OFDrfbmti2*&
   &                          ISIZE4OFDrfbmti2)
   CALL DEBUG_TGT_REAL8ARRAY('bvti1', bvti1, bvti1d, ISIZE1OFDrfbvti1&
   &                          *ISIZE2OFDrfbvti1*ISIZE3OFDrfbvti1)
   CALL DEBUG_TGT_REAL8ARRAY('bvti2', bvti2, bvti2d, ISIZE1OFDrfbvti2&
   &                          *ISIZE2OFDrfbvti2*ISIZE3OFDrfbvti2)
   CALL DEBUG_TGT_REAL8ARRAY('bmtj1', bmtj1, bmtj1d, ISIZE1OFDrfbmtj1&
   &                          *ISIZE2OFDrfbmtj1*ISIZE3OFDrfbmtj1*&
   &                          ISIZE4OFDrfbmtj1)
   CALL DEBUG_TGT_REAL8ARRAY('bmtj2', bmtj2, bmtj2d, ISIZE1OFDrfbmtj2&
   &                          *ISIZE2OFDrfbmtj2*ISIZE3OFDrfbmtj2*&
   &                          ISIZE4OFDrfbmtj2)
   CALL DEBUG_TGT_REAL8ARRAY('winf', winf, winfd, 10)
   CALL DEBUG_TGT_DISPLAY('middle')
   END IF
   bmtd => bmtk2d
   bmt => bmtk2
   bvtd => bvtk2d
   bvt => bvtk2
   END SELECT
   ! Loop over the faces of the subfaces and set the values of
   ! bmt and bvt for an implicit treatment.
   DO j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
   DO i=bcdata(nn)%icbeg,bcdata(nn)%icend
   ! Store the three components of the unit normal a bit easier.
   nnx = bcdata(nn)%norm(i, j, 1)
   nny = bcdata(nn)%norm(i, j, 2)
   nnz = bcdata(nn)%norm(i, j, 3)
   ! Determine the dot product between the outward pointing
   ! normal and the free stream velocity direction and add the
   ! possible grid velocity.
   dot = nnx*winf(ivx) + nny*winf(ivy) + nnz*winf(ivz) - bcdata(nn)%&
   &        rface(i, j)
   ! Determine whether we are dealing with an inflow or
   ! outflow boundary here.
   IF (dot .GT. zero) THEN
   ! Outflow. Simply extrapolation or zero Neumann BC
   ! of the turbulent variables.
   DO l=nt1,nt2
   bmtd(i, j, l, l) = 0.0_8
   bmt(i, j, l, l) = -one
   END DO
   ELSE
   ! Inflow. Turbulent variables are prescribed.
   DO l=nt1,nt2
   bvtd(i, j, l) = winfd(l)
   bvt(i, j, l) = winf(l)
   END DO
   END IF
   END DO
   END DO
   IF (.TRUE. .AND. DEBUG_TGT_HERE('exit', .FALSE.)) THEN
   CALL DEBUG_TGT_REAL8ARRAY('bvtj1', bvtj1, bvtj1d, ISIZE1OFDrfbvtj1*&
   &                        ISIZE2OFDrfbvtj1*ISIZE3OFDrfbvtj1)
   CALL DEBUG_TGT_REAL8ARRAY('bvtj2', bvtj2, bvtj2d, ISIZE1OFDrfbvtj2*&
   &                        ISIZE2OFDrfbvtj2*ISIZE3OFDrfbvtj2)
   CALL DEBUG_TGT_REAL8ARRAY('bmtk1', bmtk1, bmtk1d, ISIZE1OFDrfbmtk1*&
   &                        ISIZE2OFDrfbmtk1*ISIZE3OFDrfbmtk1*&
   &                        ISIZE4OFDrfbmtk1)
   CALL DEBUG_TGT_REAL8ARRAY('bmtk2', bmtk2, bmtk2d, ISIZE1OFDrfbmtk2*&
   &                        ISIZE2OFDrfbmtk2*ISIZE3OFDrfbmtk2*&
   &                        ISIZE4OFDrfbmtk2)
   CALL DEBUG_TGT_REAL8ARRAY('bvtk1', bvtk1, bvtk1d, ISIZE1OFDrfbvtk1*&
   &                        ISIZE2OFDrfbvtk1*ISIZE3OFDrfbvtk1)
   CALL DEBUG_TGT_REAL8ARRAY('bvtk2', bvtk2, bvtk2d, ISIZE1OFDrfbvtk2*&
   &                        ISIZE2OFDrfbvtk2*ISIZE3OFDrfbvtk2)
   CALL DEBUG_TGT_REAL8ARRAY('bmti1', bmti1, bmti1d, ISIZE1OFDrfbmti1*&
   &                        ISIZE2OFDrfbmti1*ISIZE3OFDrfbmti1*&
   &                        ISIZE4OFDrfbmti1)
   CALL DEBUG_TGT_REAL8ARRAY('bmti2', bmti2, bmti2d, ISIZE1OFDrfbmti2*&
   &                        ISIZE2OFDrfbmti2*ISIZE3OFDrfbmti2*&
   &                        ISIZE4OFDrfbmti2)
   CALL DEBUG_TGT_REAL8ARRAY('bvti1', bvti1, bvti1d, ISIZE1OFDrfbvti1*&
   &                        ISIZE2OFDrfbvti1*ISIZE3OFDrfbvti1)
   CALL DEBUG_TGT_REAL8ARRAY('bvti2', bvti2, bvti2d, ISIZE1OFDrfbvti2*&
   &                        ISIZE2OFDrfbvti2*ISIZE3OFDrfbvti2)
   CALL DEBUG_TGT_REAL8ARRAY('bmtj1', bmtj1, bmtj1d, ISIZE1OFDrfbmtj1*&
   &                        ISIZE2OFDrfbmtj1*ISIZE3OFDrfbmtj1*&
   &                        ISIZE4OFDrfbmtj1)
   CALL DEBUG_TGT_REAL8ARRAY('bmtj2', bmtj2, bmtj2d, ISIZE1OFDrfbmtj2*&
   &                        ISIZE2OFDrfbmtj2*ISIZE3OFDrfbmtj2*&
   &                        ISIZE4OFDrfbmtj2)
   CALL DEBUG_TGT_DISPLAY('exit')
   END IF
   END SUBROUTINE BCTURBFARFIELD_T
