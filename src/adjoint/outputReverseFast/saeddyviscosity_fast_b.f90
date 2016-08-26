!        generated by tapenade     (inria, tropics team)
!  tapenade 3.10 (r5363) -  9 sep 2014 09:53
!
!  differentiation of saeddyviscosity in reverse (adjoint) mode (with options i4 dr8 r8 noisize):
!   gradient     of useful results: *rev *w *rlv
!   with respect to varying inputs: *rev *w *rlv
!   rw status of diff variables: *rev:in-out *w:incr *rlv:incr
!   plus diff mem management of: rev:in w:in rlv:in
!      ==================================================================
!      ==================================================================
subroutine saeddyviscosity_fast_b()
!
!      ******************************************************************
!      *                                                                *
!      * saeddyviscosity computes the eddy-viscosity according to the   *
!      * spalart-allmaras model for the block given in blockpointers.   *
!      * this routine for both the original version as well as the      *
!      * modified version according to edwards.                         *
!      *                                                                *
!      ******************************************************************
!
  use constants
  use blockpointers
  use constants
  use paramturb
  implicit none
!
!      local variables.
!
  integer(kind=inttype) :: i, j, k, ii
  real(kind=realtype) :: chi, chi3, fv1, rnusa, cv13
  real(kind=realtype) :: chid, chi3d, fv1d, rnusad
  intrinsic mod
  real(kind=realtype) :: tempd
  real(kind=realtype) :: tempd0
!
!      ******************************************************************
!      *                                                                *
!      * begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
! store the cv1^3; cv1 is a constant of the spalart-allmaras model.
  cv13 = rsacv1**3
  do ii=0,ie*je*ke-1
    i = mod(ii, ie) + 1
    j = mod(ii/ie, je) + 1
    k = ii/(ie*je) + 1
    rnusa = w(i, j, k, itu1)*w(i, j, k, irho)
    chi = rnusa/rlv(i, j, k)
    chi3 = chi**3
    fv1 = chi3/(chi3+cv13)
    fv1d = rnusa*revd(i, j, k)
    tempd0 = fv1d/(cv13+chi3)
    chi3d = (1.0_8-chi3/(cv13+chi3))*tempd0
    chid = 3*chi**2*chi3d
    tempd = chid/rlv(i, j, k)
    rnusad = tempd + fv1*revd(i, j, k)
    revd(i, j, k) = 0.0_8
    rlvd(i, j, k) = rlvd(i, j, k) - rnusa*tempd/rlv(i, j, k)
    wd(i, j, k, itu1) = wd(i, j, k, itu1) + w(i, j, k, irho)*rnusad
    wd(i, j, k, irho) = wd(i, j, k, irho) + w(i, j, k, itu1)*rnusad
  end do
end subroutine saeddyviscosity_fast_b
