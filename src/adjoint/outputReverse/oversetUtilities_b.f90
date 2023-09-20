!        generated by tapenade     (inria, ecuador team)
!  tapenade 3.16 (develop) - 13 sep 2023 12:36
!
module oversetutilities_b
  implicit none

contains
!  differentiation of fractoweights in reverse (adjoint) mode (with options noisize i4 dr8 r8):
!   gradient     of useful results: weights
!   with respect to varying inputs: weights frac
!   rw status of diff variables: weights:in-out frac:out
! --------------------------------------------------
!           tapenade routine below this point
! --------------------------------------------------
  subroutine fractoweights_b(frac, fracd, weights, weightsd)
    use constants
    implicit none
    real(kind=realtype), dimension(3), intent(in) :: frac
    real(kind=realtype), dimension(3) :: fracd
    real(kind=realtype), dimension(8) :: weights
    real(kind=realtype), dimension(8) :: weightsd
    real(kind=realtype) :: tempd
    tempd = (one-frac(1))*weightsd(7)
    fracd = 0.0_8
    fracd(2) = fracd(2) + frac(1)*frac(3)*weightsd(8) + frac(3)*tempd
    fracd(3) = fracd(3) + frac(1)*frac(2)*weightsd(8) + frac(2)*tempd
    tempd = (one-frac(2))*weightsd(6)
    fracd(1) = fracd(1) + frac(2)*frac(3)*weightsd(8) + frac(3)*tempd - &
&     frac(2)*frac(3)*weightsd(7)
    weightsd(8) = 0.0_8
    weightsd(7) = 0.0_8
    fracd(3) = fracd(3) + frac(1)*tempd
    tempd = (one-frac(2))*weightsd(5)
    fracd(1) = fracd(1) - frac(3)*tempd
    fracd(3) = fracd(3) + (one-frac(1))*tempd - frac(1)*frac(2)*weightsd&
&     (4) - (one-frac(1))*frac(2)*weightsd(3) - frac(1)*(one-frac(2))*&
&     weightsd(2) - (one-frac(1))*(one-frac(2))*weightsd(1)
    tempd = (one-frac(3))*weightsd(4)
    fracd(2) = fracd(2) + frac(1)*tempd - frac(1)*frac(3)*weightsd(6) - &
&     (one-frac(1))*frac(3)*weightsd(5)
    weightsd(6) = 0.0_8
    weightsd(5) = 0.0_8
    weightsd(4) = 0.0_8
    fracd(1) = fracd(1) + frac(2)*tempd
    tempd = (one-frac(3))*weightsd(3)
    weightsd(3) = 0.0_8
    fracd(1) = fracd(1) - frac(2)*tempd
    fracd(2) = fracd(2) + (one-frac(1))*tempd
    tempd = (one-frac(3))*weightsd(2)
    weightsd(2) = 0.0_8
    fracd(1) = fracd(1) + (one-frac(2))*tempd
    fracd(2) = fracd(2) - frac(1)*tempd
    tempd = (one-frac(3))*weightsd(1)
    weightsd(1) = 0.0_8
    fracd(1) = fracd(1) - (one-frac(2))*tempd
    fracd(2) = fracd(2) - (one-frac(1))*tempd
  end subroutine fractoweights_b

! --------------------------------------------------
!           tapenade routine below this point
! --------------------------------------------------
  subroutine fractoweights(frac, weights)
    use constants
    implicit none
    real(kind=realtype), dimension(3), intent(in) :: frac
    real(kind=realtype), dimension(8), intent(out) :: weights
    weights(1) = (one-frac(1))*(one-frac(2))*(one-frac(3))
    weights(2) = frac(1)*(one-frac(2))*(one-frac(3))
    weights(3) = (one-frac(1))*frac(2)*(one-frac(3))
    weights(4) = frac(1)*frac(2)*(one-frac(3))
    weights(5) = (one-frac(1))*(one-frac(2))*frac(3)
    weights(6) = frac(1)*(one-frac(2))*frac(3)
    weights(7) = (one-frac(1))*frac(2)*frac(3)
    weights(8) = frac(1)*frac(2)*frac(3)
  end subroutine fractoweights

  subroutine fractoweights2(frac, weights)
    use constants
    implicit none
    real(kind=realtype), dimension(3), intent(in) :: frac
    real(kind=realtype), dimension(8), intent(out) :: weights
    weights(1) = (one-frac(1))*(one-frac(2))*(one-frac(3))
    weights(2) = frac(1)*(one-frac(2))*(one-frac(3))
    weights(3) = frac(1)*frac(2)*(one-frac(3))
    weights(4) = (one-frac(1))*frac(2)*(one-frac(3))
    weights(5) = (one-frac(1))*(one-frac(2))*frac(3)
    weights(6) = frac(1)*(one-frac(2))*frac(3)
    weights(7) = frac(1)*frac(2)*frac(3)
    weights(8) = (one-frac(1))*frac(2)*frac(3)
  end subroutine fractoweights2

!  differentiation of newtonupdate in reverse (adjoint) mode (with options noisize i4 dr8 r8):
!   gradient     of useful results: blk frac
!   with respect to varying inputs: xcen blk frac
!   rw status of diff variables: xcen:out blk:incr frac:in-out
  subroutine newtonupdate_b(xcen, xcend, blk, blkd, frac0, frac, fracd)
! this routine performs the newton update to recompute the new
! "frac" (u,v,w) for the point xcen. the actual search is performed
! on the the dual cell formed by the cell centers of the 3x3x3 block
! of primal nodes. this routine is ad'd with tapenade in both
! forward and reverse.
    use constants
    implicit none
! input
    real(kind=realtype), dimension(3), intent(in) :: xcen
    real(kind=realtype), dimension(3) :: xcend
    real(kind=realtype), dimension(3, 3, 3, 3), intent(in) :: blk
    real(kind=realtype), dimension(3, 3, 3, 3) :: blkd
    real(kind=realtype), dimension(3), intent(in) :: frac0
! output
    real(kind=realtype), dimension(3) :: frac
    real(kind=realtype), dimension(3) :: fracd
! working
    real(kind=realtype), dimension(3, 8) :: xn
    real(kind=realtype), dimension(3, 8) :: xnd
    real(kind=realtype) :: u, v, w, uv, uw, vw, wvu, du, dv, dw
    real(kind=realtype) :: ud, vd, wd, uvd, uwd, vwd, wvud, dud, dvd, &
&   dwd
    real(kind=realtype) :: a11, a12, a13, a21, a22, a23, a31, a32, a33, &
&   val
    real(kind=realtype) :: a11d, a12d, a13d, a21d, a22d, a23d, a31d, &
&   a32d, a33d, vald
    real(kind=realtype) :: f(3), x(3)
    real(kind=realtype) :: fd(3), xd(3)
    integer(kind=inttype), dimension(8), parameter :: indices=(/1, 2, 4&
&     , 3, 5, 6, 8, 7/)
    integer(kind=inttype) :: i, j, k, ii, ll
    real(kind=realtype), parameter :: adteps=1.e-25_realtype
    real(kind=realtype), parameter :: thresconv=1.e-10_realtype
    intrinsic sign
    intrinsic abs
    intrinsic max
    intrinsic sqrt
    real(kind=realtype) :: x1
    real(kind=realtype) :: x1d
    real(kind=realtype) :: max1
    real(kind=realtype) :: max1d
    real(kind=realtype), dimension(3) :: tempd
    real(kind=realtype), dimension(3) :: tmp
    real(kind=realtype), dimension(3) :: tmpd
    real(kind=realtype) :: temp
    real(kind=realtype) :: tempd0
    real(kind=realtype) :: temp0
    real(kind=realtype) :: tempd1
    real(kind=realtype) :: temp1
    real(kind=realtype) :: tempd2
    real(kind=realtype) :: temp2
    real(kind=realtype) :: tempd3
    integer :: branch
    integer :: ad_count
    integer :: i0
! compute the cell center locations for the 8 nodes describing the
! dual cell. note that this must be counter-clockwise ordering.
    ii = 0
    do k=1,2
      do j=1,2
        do i=1,2
          call pushinteger4(ii)
          ii = ii + 1
          xn(:, indices(ii)) = eighth*(blk(i, j, k, :)+blk(i+1, j, k, :)&
&           +blk(i, j+1, k, :)+blk(i+1, j+1, k, :)+blk(i, j, k+1, :)+blk&
&           (i+1, j, k+1, :)+blk(i, j+1, k+1, :)+blk(i+1, j+1, k+1, :))
        end do
      end do
    end do
! compute the coordinates relative to node 1.
    do i=2,8
      tmp(:) = xn(:, i) - xn(:, 1)
      xn(:, i) = tmp(:)
    end do
! compute the location of our seach point relative to the first node.
    x = xcen - xn(:, 1)
! modify the coordinates of node 3, 6, 8 and 7 such that
! they correspond to the weights of the u*v, u*w, v*w and
! u*v*w term in the transformation respectively.
    xn(1, 7) = xn(1, 7) + xn(1, 2) + xn(1, 4) + xn(1, 5) - xn(1, 3) - xn&
&     (1, 6) - xn(1, 8)
    xn(2, 7) = xn(2, 7) + xn(2, 2) + xn(2, 4) + xn(2, 5) - xn(2, 3) - xn&
&     (2, 6) - xn(2, 8)
    xn(3, 7) = xn(3, 7) + xn(3, 2) + xn(3, 4) + xn(3, 5) - xn(3, 3) - xn&
&     (3, 6) - xn(3, 8)
    xn(1, 3) = xn(1, 3) - xn(1, 2) - xn(1, 4)
    xn(2, 3) = xn(2, 3) - xn(2, 2) - xn(2, 4)
    xn(3, 3) = xn(3, 3) - xn(3, 2) - xn(3, 4)
    xn(1, 6) = xn(1, 6) - xn(1, 2) - xn(1, 5)
    xn(2, 6) = xn(2, 6) - xn(2, 2) - xn(2, 5)
    xn(3, 6) = xn(3, 6) - xn(3, 2) - xn(3, 5)
    xn(1, 8) = xn(1, 8) - xn(1, 4) - xn(1, 5)
    xn(2, 8) = xn(2, 8) - xn(2, 4) - xn(2, 5)
    xn(3, 8) = xn(3, 8) - xn(3, 4) - xn(3, 5)
! set the starting values of u, v and w based on our previous values
    u = frac0(1)
    v = frac0(2)
    w = frac0(3)
    ad_count = 1
! the newton algorithm to determine the parametric
! weights u, v and w for the given coordinate.
newtonhexa:do ll=1,15
! compute the rhs.
      uv = u*v
      uw = u*w
      vw = v*w
      wvu = u*v*w
      call pushreal8(f(1))
      f(1) = xn(1, 2)*u + xn(1, 4)*v + xn(1, 5)*w + xn(1, 3)*uv + xn(1, &
&       6)*uw + xn(1, 8)*vw + xn(1, 7)*wvu - x(1)
      call pushreal8(f(2))
      f(2) = xn(2, 2)*u + xn(2, 4)*v + xn(2, 5)*w + xn(2, 3)*uv + xn(2, &
&       6)*uw + xn(2, 8)*vw + xn(2, 7)*wvu - x(2)
      call pushreal8(f(3))
      f(3) = xn(3, 2)*u + xn(3, 4)*v + xn(3, 5)*w + xn(3, 3)*uv + xn(3, &
&       6)*uw + xn(3, 8)*vw + xn(3, 7)*wvu - x(3)
! compute the jacobian.
      a11 = xn(1, 2) + xn(1, 3)*v + xn(1, 6)*w + xn(1, 7)*vw
      a12 = xn(1, 4) + xn(1, 3)*u + xn(1, 8)*w + xn(1, 7)*uw
      a13 = xn(1, 5) + xn(1, 6)*u + xn(1, 8)*v + xn(1, 7)*uv
      a21 = xn(2, 2) + xn(2, 3)*v + xn(2, 6)*w + xn(2, 7)*vw
      a22 = xn(2, 4) + xn(2, 3)*u + xn(2, 8)*w + xn(2, 7)*uw
      a23 = xn(2, 5) + xn(2, 6)*u + xn(2, 8)*v + xn(2, 7)*uv
      a31 = xn(3, 2) + xn(3, 3)*v + xn(3, 6)*w + xn(3, 7)*vw
      a32 = xn(3, 4) + xn(3, 3)*u + xn(3, 8)*w + xn(3, 7)*uw
      a33 = xn(3, 5) + xn(3, 6)*u + xn(3, 8)*v + xn(3, 7)*uv
! compute the determinant. make sure that it is not zero
! and invert the value. the cut off is needed to be able
! to handle exceptional cases for degenerate elements.
      val = a11*(a22*a33-a32*a23) + a21*(a13*a32-a12*a33) + a31*(a12*a23&
&       -a13*a22)
      if (val .ge. 0.) then
        x1 = val
        call pushcontrol1b(0)
      else
        x1 = -val
        call pushcontrol1b(1)
      end if
      if (x1 .lt. adteps) then
        call pushreal8(max1)
        max1 = adteps
        call pushcontrol1b(0)
      else
        call pushreal8(max1)
        max1 = x1
        call pushcontrol1b(1)
      end if
      call pushreal8(val)
      val = sign(one, val)/max1
! compute the new values of u, v and w.
      du = val*((a22*a33-a23*a32)*f(1)+(a13*a32-a12*a33)*f(2)+(a12*a23-&
&       a13*a22)*f(3))
      dv = val*((a23*a31-a21*a33)*f(1)+(a11*a33-a13*a31)*f(2)+(a13*a21-&
&       a11*a23)*f(3))
      dw = val*((a21*a32-a22*a31)*f(1)+(a12*a31-a11*a32)*f(2)+(a11*a22-&
&       a12*a21)*f(3))
      call pushreal8(u)
      u = u - du
      call pushreal8(v)
      v = v - dv
      call pushreal8(w)
      w = w - dw
! exit the loop if the update of the parametric
! weights is below the threshold
      call pushreal8(val)
      val = sqrt(du*du + dv*dv + dw*dw)
      if (val .le. thresconv) then
        goto 100
      else
        ad_count = ad_count + 1
      end if
    end do newtonhexa
    call pushcontrol1b(0)
    call pushinteger4(ad_count)
    goto 110
 100 call pushcontrol1b(1)
    call pushinteger4(ad_count)
 110 wd = fracd(3)
    fracd(3) = 0.0_8
    vd = fracd(2)
    fracd(2) = 0.0_8
    ud = fracd(1)
    fracd(1) = 0.0_8
    call popinteger4(ad_count)
    do 120 i0=1,ad_count
      if (i0 .eq. 1) then
        call popcontrol1b(branch)
        if (branch .eq. 0) then
          fd = 0.0_8
          xd = 0.0_8
          xnd = 0.0_8
          goto 120
        else
          fd = 0.0_8
          xd = 0.0_8
          xnd = 0.0_8
        end if
      end if
      call popreal8(val)
      call popreal8(w)
      dwd = -wd
      call popreal8(v)
      dvd = -vd
      call popreal8(u)
      dud = -ud
      vw = v*w
      a21 = xn(2, 2) + xn(2, 3)*v + xn(2, 6)*w + xn(2, 7)*vw
      uw = u*w
      a22 = xn(2, 4) + xn(2, 3)*u + xn(2, 8)*w + xn(2, 7)*uw
      a31 = xn(3, 2) + xn(3, 3)*v + xn(3, 6)*w + xn(3, 7)*vw
      a32 = xn(3, 4) + xn(3, 3)*u + xn(3, 8)*w + xn(3, 7)*uw
      a11 = xn(1, 2) + xn(1, 3)*v + xn(1, 6)*w + xn(1, 7)*vw
      a12 = xn(1, 4) + xn(1, 3)*u + xn(1, 8)*w + xn(1, 7)*uw
      temp2 = a11*a22 - a12*a21
      temp1 = a12*a31 - a11*a32
      temp0 = a21*a32 - a22*a31
      vald = (temp0*f(1)+temp1*f(2)+temp2*f(3))*dwd
      tempd0 = val*dwd
      tempd1 = f(1)*tempd0
      fd(1) = fd(1) + temp0*tempd0
      tempd2 = f(2)*tempd0
      fd(2) = fd(2) + temp1*tempd0
      tempd3 = f(3)*tempd0
      fd(3) = fd(3) + temp2*tempd0
      a11d = a22*tempd3 - a32*tempd2
      a22d = a11*tempd3 - a31*tempd1
      a12d = a31*tempd2 - a21*tempd3
      a21d = a32*tempd1 - a12*tempd3
      a31d = a12*tempd2 - a22*tempd1
      a32d = a21*tempd1 - a11*tempd2
      uv = u*v
      a23 = xn(2, 5) + xn(2, 6)*u + xn(2, 8)*v + xn(2, 7)*uv
      a33 = xn(3, 5) + xn(3, 6)*u + xn(3, 8)*v + xn(3, 7)*uv
      a13 = xn(1, 5) + xn(1, 6)*u + xn(1, 8)*v + xn(1, 7)*uv
      temp2 = a13*a21 - a11*a23
      temp1 = a11*a33 - a13*a31
      temp0 = a23*a31 - a21*a33
      vald = vald + (temp0*f(1)+temp1*f(2)+temp2*f(3))*dvd
      tempd0 = val*dvd
      tempd1 = f(1)*tempd0
      tempd2 = f(2)*tempd0
      tempd3 = f(3)*tempd0
      a13d = a21*tempd3 - a31*tempd2
      a21d = a21d + a13*tempd3 - a33*tempd1
      a11d = a11d + a33*tempd2 - a23*tempd3
      a23d = -(a11*tempd3)
      a33d = a11*tempd2 - a21*tempd1
      a31d = a31d + a23*tempd1 - a13*tempd2
      temp = a22*a33 - a23*a32
      tempd3 = val*dud
      fd(1) = fd(1) + temp0*tempd0 + temp*tempd3
      temp0 = a13*a32 - a12*a33
      fd(2) = fd(2) + temp1*tempd0 + temp0*tempd3
      temp1 = a12*a23 - a13*a22
      fd(3) = fd(3) + temp2*tempd0 + temp1*tempd3
      vald = vald + (temp*f(1)+temp0*f(2)+temp1*f(3))*dud
      tempd0 = f(1)*tempd3
      tempd2 = f(3)*tempd3
      a23d = a23d + a31*tempd1 + a12*tempd2 - a32*tempd0
      tempd1 = f(2)*tempd3
      a12d = a12d + a23*tempd2 - a33*tempd1
      a13d = a13d + a32*tempd1 - a22*tempd2
      a22d = a22d + a33*tempd0 - a13*tempd2
      a32d = a32d + a13*tempd1 - a23*tempd0
      a33d = a33d + a22*tempd0 - a12*tempd1
      call popreal8(val)
      max1d = -(sign(one, val)*vald/max1**2)
      call popcontrol1b(branch)
      if (branch .eq. 0) then
        call popreal8(max1)
        x1d = 0.0_8
      else
        call popreal8(max1)
        x1d = max1d
      end if
      call popcontrol1b(branch)
      if (branch .eq. 0) then
        vald = x1d
      else
        vald = -x1d
      end if
      wvud = xn(3, 7)*fd(3) + xn(2, 7)*fd(2) + xn(1, 7)*fd(1)
      wvu = u*v*w
      a11d = a11d + (a22*a33-a32*a23)*vald
      tempd0 = a11*vald
      a21d = a21d + (a13*a32-a12*a33)*vald
      tempd1 = a21*vald
      a31d = a31d + (a12*a23-a13*a22)*vald
      tempd2 = a31*vald
      a12d = a12d + a23*tempd2 - a33*tempd1
      a23d = a23d + a12*tempd2 - a32*tempd0
      a13d = a13d + a32*tempd1 - a22*tempd2
      a22d = a22d + a33*tempd0 - a13*tempd2
      a32d = a32d + a13*tempd1 - a23*tempd0
      a33d = a33d + a22*tempd0 - a12*tempd1
      xnd(3, 5) = xnd(3, 5) + a33d + w*fd(3)
      xnd(3, 6) = xnd(3, 6) + u*a33d + w*a31d + uw*fd(3)
      xnd(3, 8) = xnd(3, 8) + v*a33d + w*a32d + vw*fd(3)
      xnd(3, 7) = xnd(3, 7) + uv*a33d + uw*a32d + vw*a31d + wvu*fd(3)
      uvd = xn(3, 7)*a33d + xn(2, 7)*a23d + xn(1, 7)*a13d + xn(3, 3)*fd(&
&       3) + xn(2, 3)*fd(2) + xn(1, 3)*fd(1)
      xnd(3, 4) = xnd(3, 4) + a32d + v*fd(3)
      xnd(3, 3) = xnd(3, 3) + u*a32d + v*a31d + uv*fd(3)
      uwd = xn(3, 7)*a32d + xn(2, 7)*a22d + xn(1, 7)*a12d + xn(3, 6)*fd(&
&       3) + xn(2, 6)*fd(2) + xn(1, 6)*fd(1)
      ud = ud + xn(3, 6)*a33d + xn(3, 3)*a32d + xn(2, 6)*a23d + xn(2, 3)&
&       *a22d + xn(1, 6)*a13d + xn(1, 3)*a12d + xn(3, 2)*fd(3) + xn(2, 2&
&       )*fd(2) + xn(1, 2)*fd(1) + v*w*wvud + w*uwd + v*uvd
      xnd(3, 2) = xnd(3, 2) + a31d + u*fd(3)
      vwd = xn(3, 7)*a31d + xn(2, 7)*a21d + xn(1, 7)*a11d + xn(3, 8)*fd(&
&       3) + xn(2, 8)*fd(2) + xn(1, 8)*fd(1)
      vd = vd + xn(3, 8)*a33d + xn(3, 3)*a31d + xn(2, 8)*a23d + xn(2, 3)&
&       *a21d + xn(1, 8)*a13d + xn(1, 3)*a11d + xn(3, 4)*fd(3) + xn(2, 4&
&       )*fd(2) + xn(1, 4)*fd(1) + u*w*wvud + w*vwd + u*uvd
      wd = wd + xn(3, 8)*a32d + xn(3, 6)*a31d + xn(2, 8)*a22d + xn(2, 6)&
&       *a21d + xn(1, 8)*a12d + xn(1, 6)*a11d + xn(3, 5)*fd(3) + xn(2, 5&
&       )*fd(2) + xn(1, 5)*fd(1) + u*v*wvud + v*vwd + u*uwd
      xnd(2, 5) = xnd(2, 5) + a23d + w*fd(2)
      xnd(2, 6) = xnd(2, 6) + u*a23d + w*a21d + uw*fd(2)
      xnd(2, 8) = xnd(2, 8) + v*a23d + w*a22d + vw*fd(2)
      xnd(2, 7) = xnd(2, 7) + uv*a23d + uw*a22d + vw*a21d + wvu*fd(2)
      xnd(2, 4) = xnd(2, 4) + a22d + v*fd(2)
      xnd(2, 3) = xnd(2, 3) + u*a22d + v*a21d + uv*fd(2)
      xnd(2, 2) = xnd(2, 2) + a21d + u*fd(2)
      xnd(1, 5) = xnd(1, 5) + a13d + w*fd(1)
      xnd(1, 6) = xnd(1, 6) + u*a13d + w*a11d + uw*fd(1)
      xnd(1, 8) = xnd(1, 8) + v*a13d + w*a12d + vw*fd(1)
      xnd(1, 7) = xnd(1, 7) + uv*a13d + uw*a12d + vw*a11d + wvu*fd(1)
      xnd(1, 4) = xnd(1, 4) + a12d + v*fd(1)
      xnd(1, 3) = xnd(1, 3) + u*a12d + v*a11d + uv*fd(1)
      xnd(1, 2) = xnd(1, 2) + a11d + u*fd(1)
      call popreal8(f(3))
      xd(3) = xd(3) - fd(3)
      fd(3) = 0.0_8
      call popreal8(f(2))
      xd(2) = xd(2) - fd(2)
      fd(2) = 0.0_8
      call popreal8(f(1))
      xd(1) = xd(1) - fd(1)
      fd(1) = 0.0_8
 120 continue
    xnd(3, 4) = xnd(3, 4) + xnd(3, 7) - xnd(3, 8) - xnd(3, 3)
    xnd(3, 5) = xnd(3, 5) + xnd(3, 7) - xnd(3, 8) - xnd(3, 6)
    xnd(2, 4) = xnd(2, 4) + xnd(2, 7) - xnd(2, 8) - xnd(2, 3)
    xnd(2, 5) = xnd(2, 5) + xnd(2, 7) - xnd(2, 8) - xnd(2, 6)
    xnd(1, 4) = xnd(1, 4) + xnd(1, 7) - xnd(1, 8) - xnd(1, 3)
    xnd(1, 5) = xnd(1, 5) + xnd(1, 7) - xnd(1, 8) - xnd(1, 6)
    xnd(3, 2) = xnd(3, 2) + xnd(3, 7) - xnd(3, 6) - xnd(3, 3)
    xnd(2, 2) = xnd(2, 2) + xnd(2, 7) - xnd(2, 6) - xnd(2, 3)
    xnd(1, 2) = xnd(1, 2) + xnd(1, 7) - xnd(1, 6) - xnd(1, 3)
    xnd(3, 3) = xnd(3, 3) - xnd(3, 7)
    xnd(3, 6) = xnd(3, 6) - xnd(3, 7)
    xnd(3, 8) = xnd(3, 8) - xnd(3, 7)
    xnd(2, 3) = xnd(2, 3) - xnd(2, 7)
    xnd(2, 6) = xnd(2, 6) - xnd(2, 7)
    xnd(2, 8) = xnd(2, 8) - xnd(2, 7)
    xnd(1, 3) = xnd(1, 3) - xnd(1, 7)
    xnd(1, 6) = xnd(1, 6) - xnd(1, 7)
    xnd(1, 8) = xnd(1, 8) - xnd(1, 7)
    xcend = 0.0_8
    xcend = xd
    xnd(:, 1) = xnd(:, 1) - xd
    do i=8,2,-1
      tmpd = xnd(:, i)
      xnd(:, i) = tmpd
      xnd(:, 1) = xnd(:, 1) - tmpd
    end do
    do k=2,1,-1
      do j=2,1,-1
        do i=2,1,-1
          tempd = eighth*xnd(:, indices(ii))
          xnd(:, indices(ii)) = 0.0_8
          blkd(i, j, k, :) = blkd(i, j, k, :) + tempd
          blkd(i+1, j, k, :) = blkd(i+1, j, k, :) + tempd
          blkd(i, j+1, k, :) = blkd(i, j+1, k, :) + tempd
          blkd(i+1, j+1, k, :) = blkd(i+1, j+1, k, :) + tempd
          blkd(i, j, k+1, :) = blkd(i, j, k+1, :) + tempd
          blkd(i+1, j, k+1, :) = blkd(i+1, j, k+1, :) + tempd
          blkd(i, j+1, k+1, :) = blkd(i, j+1, k+1, :) + tempd
          blkd(i+1, j+1, k+1, :) = blkd(i+1, j+1, k+1, :) + tempd
          call popinteger4(ii)
        end do
      end do
    end do
  end subroutine newtonupdate_b

  subroutine newtonupdate(xcen, blk, frac0, frac)
! this routine performs the newton update to recompute the new
! "frac" (u,v,w) for the point xcen. the actual search is performed
! on the the dual cell formed by the cell centers of the 3x3x3 block
! of primal nodes. this routine is ad'd with tapenade in both
! forward and reverse.
    use constants
    implicit none
! input
    real(kind=realtype), dimension(3), intent(in) :: xcen
    real(kind=realtype), dimension(3, 3, 3, 3), intent(in) :: blk
    real(kind=realtype), dimension(3), intent(in) :: frac0
! output
    real(kind=realtype), dimension(3), intent(out) :: frac
! working
    real(kind=realtype), dimension(3, 8) :: xn
    real(kind=realtype) :: u, v, w, uv, uw, vw, wvu, du, dv, dw
    real(kind=realtype) :: a11, a12, a13, a21, a22, a23, a31, a32, a33, &
&   val
    real(kind=realtype) :: f(3), x(3)
    integer(kind=inttype), dimension(8), parameter :: indices=(/1, 2, 4&
&     , 3, 5, 6, 8, 7/)
    integer(kind=inttype) :: i, j, k, ii, ll
    real(kind=realtype), parameter :: adteps=1.e-25_realtype
    real(kind=realtype), parameter :: thresconv=1.e-10_realtype
    intrinsic sign
    intrinsic abs
    intrinsic max
    intrinsic sqrt
    real(kind=realtype) :: x1
    real(kind=realtype) :: max1
! compute the cell center locations for the 8 nodes describing the
! dual cell. note that this must be counter-clockwise ordering.
    ii = 0
    do k=1,2
      do j=1,2
        do i=1,2
          ii = ii + 1
          xn(:, indices(ii)) = eighth*(blk(i, j, k, :)+blk(i+1, j, k, :)&
&           +blk(i, j+1, k, :)+blk(i+1, j+1, k, :)+blk(i, j, k+1, :)+blk&
&           (i+1, j, k+1, :)+blk(i, j+1, k+1, :)+blk(i+1, j+1, k+1, :))
        end do
      end do
    end do
! compute the coordinates relative to node 1.
    do i=2,8
      xn(:, i) = xn(:, i) - xn(:, 1)
    end do
! compute the location of our seach point relative to the first node.
    x = xcen - xn(:, 1)
! modify the coordinates of node 3, 6, 8 and 7 such that
! they correspond to the weights of the u*v, u*w, v*w and
! u*v*w term in the transformation respectively.
    xn(1, 7) = xn(1, 7) + xn(1, 2) + xn(1, 4) + xn(1, 5) - xn(1, 3) - xn&
&     (1, 6) - xn(1, 8)
    xn(2, 7) = xn(2, 7) + xn(2, 2) + xn(2, 4) + xn(2, 5) - xn(2, 3) - xn&
&     (2, 6) - xn(2, 8)
    xn(3, 7) = xn(3, 7) + xn(3, 2) + xn(3, 4) + xn(3, 5) - xn(3, 3) - xn&
&     (3, 6) - xn(3, 8)
    xn(1, 3) = xn(1, 3) - xn(1, 2) - xn(1, 4)
    xn(2, 3) = xn(2, 3) - xn(2, 2) - xn(2, 4)
    xn(3, 3) = xn(3, 3) - xn(3, 2) - xn(3, 4)
    xn(1, 6) = xn(1, 6) - xn(1, 2) - xn(1, 5)
    xn(2, 6) = xn(2, 6) - xn(2, 2) - xn(2, 5)
    xn(3, 6) = xn(3, 6) - xn(3, 2) - xn(3, 5)
    xn(1, 8) = xn(1, 8) - xn(1, 4) - xn(1, 5)
    xn(2, 8) = xn(2, 8) - xn(2, 4) - xn(2, 5)
    xn(3, 8) = xn(3, 8) - xn(3, 4) - xn(3, 5)
! set the starting values of u, v and w based on our previous values
    u = frac0(1)
    v = frac0(2)
    w = frac0(3)
! the newton algorithm to determine the parametric
! weights u, v and w for the given coordinate.
newtonhexa:do ll=1,15
! compute the rhs.
      uv = u*v
      uw = u*w
      vw = v*w
      wvu = u*v*w
      f(1) = xn(1, 2)*u + xn(1, 4)*v + xn(1, 5)*w + xn(1, 3)*uv + xn(1, &
&       6)*uw + xn(1, 8)*vw + xn(1, 7)*wvu - x(1)
      f(2) = xn(2, 2)*u + xn(2, 4)*v + xn(2, 5)*w + xn(2, 3)*uv + xn(2, &
&       6)*uw + xn(2, 8)*vw + xn(2, 7)*wvu - x(2)
      f(3) = xn(3, 2)*u + xn(3, 4)*v + xn(3, 5)*w + xn(3, 3)*uv + xn(3, &
&       6)*uw + xn(3, 8)*vw + xn(3, 7)*wvu - x(3)
! compute the jacobian.
      a11 = xn(1, 2) + xn(1, 3)*v + xn(1, 6)*w + xn(1, 7)*vw
      a12 = xn(1, 4) + xn(1, 3)*u + xn(1, 8)*w + xn(1, 7)*uw
      a13 = xn(1, 5) + xn(1, 6)*u + xn(1, 8)*v + xn(1, 7)*uv
      a21 = xn(2, 2) + xn(2, 3)*v + xn(2, 6)*w + xn(2, 7)*vw
      a22 = xn(2, 4) + xn(2, 3)*u + xn(2, 8)*w + xn(2, 7)*uw
      a23 = xn(2, 5) + xn(2, 6)*u + xn(2, 8)*v + xn(2, 7)*uv
      a31 = xn(3, 2) + xn(3, 3)*v + xn(3, 6)*w + xn(3, 7)*vw
      a32 = xn(3, 4) + xn(3, 3)*u + xn(3, 8)*w + xn(3, 7)*uw
      a33 = xn(3, 5) + xn(3, 6)*u + xn(3, 8)*v + xn(3, 7)*uv
! compute the determinant. make sure that it is not zero
! and invert the value. the cut off is needed to be able
! to handle exceptional cases for degenerate elements.
      val = a11*(a22*a33-a32*a23) + a21*(a13*a32-a12*a33) + a31*(a12*a23&
&       -a13*a22)
      if (val .ge. 0.) then
        x1 = val
      else
        x1 = -val
      end if
      if (x1 .lt. adteps) then
        max1 = adteps
      else
        max1 = x1
      end if
      val = sign(one, val)/max1
! compute the new values of u, v and w.
      du = val*((a22*a33-a23*a32)*f(1)+(a13*a32-a12*a33)*f(2)+(a12*a23-&
&       a13*a22)*f(3))
      dv = val*((a23*a31-a21*a33)*f(1)+(a11*a33-a13*a31)*f(2)+(a13*a21-&
&       a11*a23)*f(3))
      dw = val*((a21*a32-a22*a31)*f(1)+(a12*a31-a11*a32)*f(2)+(a11*a22-&
&       a12*a21)*f(3))
      u = u - du
      v = v - dv
      w = w - dw
! exit the loop if the update of the parametric
! weights is below the threshold
      val = sqrt(du*du + dv*dv + dw*dw)
      if (val .le. thresconv) goto 100
    end do newtonhexa
! we would *like* that all solutions fall inside the hexa, but we
! can't be picky here since we are not changing the donors. so
! whatever the u,v,w is we have to accept. even if it is greater than
! 1 or less than zero, it shouldn't be by much.
 100 frac(1) = u
    frac(2) = v
    frac(3) = w
  end subroutine newtonupdate

end module oversetutilities_b

