!        generated by tapenade     (inria, tropics team)
!  tapenade 3.10 (r5363) -  9 sep 2014 09:53
!
!  differentiation of metric_block in forward (tangent) mode (with options i4 dr8 r8):
!   variations   of useful results: *si *sj *sk
!   with respect to varying inputs: *x
!   plus diff mem management of: x:in si:in sj:in sk:in
subroutine metric_block_d()
  use constants
  use blockpointers
  implicit none
! local variables.
  integer(kind=inttype) :: i, j, k, n, m, l, ii
  real(kind=realtype) :: fact
  real(kind=realtype) :: xxp, yyp, zzp
  real(kind=realtype), dimension(3) :: v1, v2
  real(kind=realtype), dimension(3) :: v1d, v2d
  intrinsic mod
! set the factor in the surface normals computation. for a
! left handed block this factor is negative, such that the
! normals still point in the direction of increasing index.
! the formulae used later on assume a right handed block
! and fact is used to correct this for a left handed block,
! as well as the scaling factor of 0.5
  if (righthanded) then
    fact = half
    sid = 0.0_8
    v1d = 0.0_8
    v2d = 0.0_8
  else
    fact = -half
    sid = 0.0_8
    v1d = 0.0_8
    v2d = 0.0_8
  end if
!
! **************************************************************
! *                                                            *
! * computation of the face normals in i-, j- and k-direction. *
! * formula's are valid for a right handed block; for a left   *
! * handed block the correct orientation is obtained via fact. *
! * the normals point in the direction of increasing index.    *
! * the absolute value of fact is 0.5, because the cross       *
! * product of the two diagonals is twice the normal vector.   *
! *                                                            *
! * note that also the normals of the first level halo cells   *
! * are computed. these are needed for the viscous fluxes.     *
! *                                                            *
! **************************************************************
!
! projected areas of cell faces in the i direction.
  do ii=0,ke*je*(ie+1)-1
! 0:ie
    i = mod(ii, ie + 1) + 0
!1:je
    j = mod(ii/(ie+1), je) + 1
!1:ke
    k = ii/((ie+1)*je) + 1
    n = k - 1
    m = j - 1
! determine the two diagonal vectors of the face.
    v1d(1) = xd(i, j, n, 1) - xd(i, m, k, 1)
    v1(1) = x(i, j, n, 1) - x(i, m, k, 1)
    v1d(2) = xd(i, j, n, 2) - xd(i, m, k, 2)
    v1(2) = x(i, j, n, 2) - x(i, m, k, 2)
    v1d(3) = xd(i, j, n, 3) - xd(i, m, k, 3)
    v1(3) = x(i, j, n, 3) - x(i, m, k, 3)
    v2d(1) = xd(i, j, k, 1) - xd(i, m, n, 1)
    v2(1) = x(i, j, k, 1) - x(i, m, n, 1)
    v2d(2) = xd(i, j, k, 2) - xd(i, m, n, 2)
    v2(2) = x(i, j, k, 2) - x(i, m, n, 2)
    v2d(3) = xd(i, j, k, 3) - xd(i, m, n, 3)
    v2(3) = x(i, j, k, 3) - x(i, m, n, 3)
! the face normal, which is the cross product of the two
! diagonal vectors times fact; remember that fact is
! either -0.5 or 0.5.
    sid(i, j, k, 1) = fact*(v1d(2)*v2(3)+v1(2)*v2d(3)-v1d(3)*v2(2)-v1(3)&
&     *v2d(2))
    si(i, j, k, 1) = fact*(v1(2)*v2(3)-v1(3)*v2(2))
    sid(i, j, k, 2) = fact*(v1d(3)*v2(1)+v1(3)*v2d(1)-v1d(1)*v2(3)-v1(1)&
&     *v2d(3))
    si(i, j, k, 2) = fact*(v1(3)*v2(1)-v1(1)*v2(3))
    sid(i, j, k, 3) = fact*(v1d(1)*v2(2)+v1(1)*v2d(2)-v1d(2)*v2(1)-v1(2)&
&     *v2d(1))
    si(i, j, k, 3) = fact*(v1(1)*v2(2)-v1(2)*v2(1))
  end do
  sjd = 0.0_8
! projected areas of cell faces in the j direction
  do ii=0,ke*(je+1)*ie-1
! 1:ie
    i = mod(ii, ie) + 1
!0:je
    j = mod(ii/ie, je + 1) + 0
!1:ke
    k = ii/(ie*(je+1)) + 1
    n = k - 1
    l = i - 1
! determine the two diagonal vectors of the face.
    v1d(1) = xd(i, j, n, 1) - xd(l, j, k, 1)
    v1(1) = x(i, j, n, 1) - x(l, j, k, 1)
    v1d(2) = xd(i, j, n, 2) - xd(l, j, k, 2)
    v1(2) = x(i, j, n, 2) - x(l, j, k, 2)
    v1d(3) = xd(i, j, n, 3) - xd(l, j, k, 3)
    v1(3) = x(i, j, n, 3) - x(l, j, k, 3)
    v2d(1) = xd(l, j, n, 1) - xd(i, j, k, 1)
    v2(1) = x(l, j, n, 1) - x(i, j, k, 1)
    v2d(2) = xd(l, j, n, 2) - xd(i, j, k, 2)
    v2(2) = x(l, j, n, 2) - x(i, j, k, 2)
    v2d(3) = xd(l, j, n, 3) - xd(i, j, k, 3)
    v2(3) = x(l, j, n, 3) - x(i, j, k, 3)
! the face normal, which is the cross product of the two
! diagonal vectors times fact; remember that fact is
! either -0.5 or 0.5.
    sjd(i, j, k, 1) = fact*(v1d(2)*v2(3)+v1(2)*v2d(3)-v1d(3)*v2(2)-v1(3)&
&     *v2d(2))
    sj(i, j, k, 1) = fact*(v1(2)*v2(3)-v1(3)*v2(2))
    sjd(i, j, k, 2) = fact*(v1d(3)*v2(1)+v1(3)*v2d(1)-v1d(1)*v2(3)-v1(1)&
&     *v2d(3))
    sj(i, j, k, 2) = fact*(v1(3)*v2(1)-v1(1)*v2(3))
    sjd(i, j, k, 3) = fact*(v1d(1)*v2(2)+v1(1)*v2d(2)-v1d(2)*v2(1)-v1(2)&
&     *v2d(1))
    sj(i, j, k, 3) = fact*(v1(1)*v2(2)-v1(2)*v2(1))
  end do
  skd = 0.0_8
! projected areas of cell faces in the k direction.
  do ii=0,(ke+1)*je*ie-1
! 1:ie
    i = mod(ii, ie) + 1
!1:je
    j = mod(ii/ie, je) + 1
!0:ke
    k = ii/(ie*je) + 0
    m = j - 1
    l = i - 1
! determine the two diagonal vectors of the face.
    v1d(1) = xd(i, j, k, 1) - xd(l, m, k, 1)
    v1(1) = x(i, j, k, 1) - x(l, m, k, 1)
    v1d(2) = xd(i, j, k, 2) - xd(l, m, k, 2)
    v1(2) = x(i, j, k, 2) - x(l, m, k, 2)
    v1d(3) = xd(i, j, k, 3) - xd(l, m, k, 3)
    v1(3) = x(i, j, k, 3) - x(l, m, k, 3)
    v2d(1) = xd(l, j, k, 1) - xd(i, m, k, 1)
    v2(1) = x(l, j, k, 1) - x(i, m, k, 1)
    v2d(2) = xd(l, j, k, 2) - xd(i, m, k, 2)
    v2(2) = x(l, j, k, 2) - x(i, m, k, 2)
    v2d(3) = xd(l, j, k, 3) - xd(i, m, k, 3)
    v2(3) = x(l, j, k, 3) - x(i, m, k, 3)
! the face normal, which is the cross product of the two
! diagonal vectors times fact; remember that fact is
! either -0.5 or 0.5.
    skd(i, j, k, 1) = fact*(v1d(2)*v2(3)+v1(2)*v2d(3)-v1d(3)*v2(2)-v1(3)&
&     *v2d(2))
    sk(i, j, k, 1) = fact*(v1(2)*v2(3)-v1(3)*v2(2))
    skd(i, j, k, 2) = fact*(v1d(3)*v2(1)+v1(3)*v2d(1)-v1d(1)*v2(3)-v1(1)&
&     *v2d(3))
    sk(i, j, k, 2) = fact*(v1(3)*v2(1)-v1(1)*v2(3))
    skd(i, j, k, 3) = fact*(v1d(1)*v2(2)+v1(1)*v2d(2)-v1d(2)*v2(1)-v1(2)&
&     *v2d(1))
    sk(i, j, k, 3) = fact*(v1(1)*v2(2)-v1(2)*v2(1))
  end do
end subroutine metric_block_d
