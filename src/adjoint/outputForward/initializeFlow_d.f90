!        generated by tapenade     (inria, ecuador team)
!  tapenade 3.16 (develop) - 10 nov 2023 18:24
!
module initializeflow_d
  use constants, only : inttype, realtype, maxstringlen
  implicit none
  save 

contains
!  differentiation of referencestate in forward (tangent) mode (with options i4 dr8 r8):
!   variations   of useful results: pinf timeref rhoinf muref tref
!                winf muinf uinf pinfcorr pref rgas rhoref
!   with respect to varying inputs: tinfdim rhoinfdim pinfdim mach
!                veldirfreestream machcoef
!   rw status of diff variables: tinfdim:in pinf:out timeref:out
!                rhoinf:out muref:out rhoinfdim:in tref:out winf:out
!                muinf:out uinf:out pinfcorr:out pref:out pinfdim:in
!                rgas:out rhoref:out mach:in veldirfreestream:in
!                machcoef:in
  subroutine referencestate_d()
!
!       the original version has been nuked since the computations are
!       no longer necessary when calling from python
!       this is the most compliclated routine in all of adflow. it is
!       stupidly complicated. this is most likely the reason your
!       derivatives are wrong. you don't understand this routine
!       and its effects.
!       this routine *requries* the following as input:
!       mach, pinfdim, tinfdim, rhoinfdim, rgasdim (machcoef non-sa
!        turbulence only)
!       optionally, pref, rhoref and tref are used if they are
!       are non-negative. this only happens when you want the equations
!       normalized by values other than the freestream
!      * this routine computes as output:
!      *   muinfdim, (unused anywhere in code)
!         pref, rhoref, tref, muref, timeref ('dimensional' reference)
!         pinf, pinfcorr, rhoinf, uinf, rgas, muinf, gammainf and winf
!         (non-dimensionalized values used in actual computations)
!
    use constants
    use paramturb
    use inputphysics, only : equations, mach, machd, machcoef, &
&   machcoefd, musuthdim, tsuthdim, veldirfreestream, veldirfreestreamd,&
&   rgasdim, ssuthdim, eddyvisinfratio, turbmodel, turbintensityinf
    use flowvarrefstate, only : pinfdim, pinfdimd, tinfdim, tinfdimd, &
&   rhoinfdim, rhoinfdimd, muinfdim, muinfdimd, pref, prefd, rhoref, &
&   rhorefd, tref, trefd, muref, murefd, timeref, timerefd, uref, urefd,&
&   href, hrefd, pinf, pinfd, pinfcorr, pinfcorrd, rhoinf, rhoinfd, uinf&
&   , uinfd, rgas, rgasd, muinf, muinfd, gammainf, winf, winfd, nw, nwf,&
&   kpresent, winf, winfd
    use flowutils_d, only : computegamma, etot, etot_d
    use turbutils_d, only : sanuknowneddyratio, sanuknowneddyratio_d
    implicit none
    integer(kind=inttype) :: sps, nn, mm, ierr
    real(kind=realtype) :: gm1, ratio
    real(kind=realtype) :: nuinf, ktmp, uinf2
    real(kind=realtype) :: nuinfd, ktmpd, uinf2d
    real(kind=realtype) :: vinf, zinf, tmp1(1), tmp2(1)
    real(kind=realtype) :: vinfd, zinfd
    intrinsic sqrt
    real(kind=realtype) :: arg1
    real(kind=realtype) :: arg1d
    real(kind=realtype) :: result1
    real(kind=realtype) :: result1d
    real(kind=realtype) :: temp
    real(kind=realtype) :: temp0
! compute the dimensional viscosity from sutherland's law
    temp = (tinfdim/tsuthdim)**1.5_realtype/(ssuthdim+tinfdim)
    muinfdimd = musuthdim*(tsuthdim+ssuthdim)*(1.5_realtype*(tinfdim/&
&     tsuthdim)**0.5/tsuthdim-temp)*tinfdimd/(ssuthdim+tinfdim)
    muinfdim = musuthdim*(tsuthdim+ssuthdim)*temp
! set the reference values. they *could* be different from the
! free-stream values for an internal flow simulation. for now,
! we just use the actual free stream values.
    prefd = pinfdimd
    pref = pinfdim
    trefd = tinfdimd
    tref = tinfdim
    rhorefd = rhoinfdimd
    rhoref = rhoinfdim
! compute the value of muref, such that the nondimensional
! equations are identical to the dimensional ones.
! note that in the non-dimensionalization of muref there is
! a reference length. however this reference length is 1.0
! in this code, because the coordinates are converted to
! meters.
    temp = sqrt(pref*rhoref)
    if (pref*rhoref .eq. 0.0_8) then
      murefd = 0.0_8
    else
      murefd = (rhoref*prefd+pref*rhorefd)/(2.0*temp)
    end if
    muref = temp
! compute timeref for a correct nondimensionalization of the
! unsteady equations. some story as for the reference viscosity
! concerning the reference length.
    temp = rhoref/pref
    temp0 = sqrt(temp)
    if (temp .eq. 0.0_8) then
      timerefd = 0.0_8
    else
      timerefd = (rhorefd-temp*prefd)/(2.0*temp0*pref)
    end if
    timeref = temp0
    href = pref/rhoref
    uref = sqrt(href)
! compute the nondimensional pressure, density, velocity,
! viscosity and gas constant.
    pinfd = (pinfdimd-pinfdim*prefd/pref)/pref
    pinf = pinfdim/pref
    rhoinfd = (rhoinfdimd-rhoinfdim*rhorefd/rhoref)/rhoref
    rhoinf = rhoinfdim/rhoref
    arg1d = gammainf*(pinfd-pinf*rhoinfd/rhoinf)/rhoinf
    arg1 = gammainf*pinf/rhoinf
    temp0 = sqrt(arg1)
    if (arg1 .eq. 0.0_8) then
      result1d = 0.0_8
    else
      result1d = arg1d/(2.0*temp0)
    end if
    result1 = temp0
    uinfd = result1*machd + mach*result1d
    uinf = mach*result1
    temp0 = rhoref*tref/pref
    rgasd = rgasdim*(tref*rhorefd+rhoref*trefd-temp0*prefd)/pref
    rgas = rgasdim*temp0
    muinfd = (muinfdimd-muinfdim*murefd/muref)/muref
    muinf = muinfdim/muref
    tmp1(1) = tinfdim
    call computegamma(tmp1, tmp2, 1)
    gammainf = tmp2(1)
! ----------------------------------------
!      compute the final winf
! ----------------------------------------
! allocate the memory for winf if necessary
! zero out the winf first
    winf(:) = zero
! set the reference value of the flow variables, except the total
! energy. this will be computed at the end of this routine.
    winfd = 0.0_8
    winfd(irho) = rhoinfd
    winf(irho) = rhoinf
    winfd(ivx) = veldirfreestream(1)*uinfd + uinf*veldirfreestreamd(1)
    winf(ivx) = uinf*veldirfreestream(1)
    winfd(ivy) = veldirfreestream(2)*uinfd + uinf*veldirfreestreamd(2)
    winf(ivy) = uinf*veldirfreestream(2)
    winfd(ivz) = veldirfreestream(3)*uinfd + uinf*veldirfreestreamd(3)
    winf(ivz) = uinf*veldirfreestream(3)
! compute the velocity squared based on machcoef. this gives a
! better indication of the 'speed' of the flow so the turubulence
! intensity ration is more meaningful especially for moving
! geometries. (not used in sa model)
    temp0 = pinf/rhoinf
    uinf2d = gammainf*(temp0*2*machcoef*machcoefd+machcoef**2*(pinfd-&
&     temp0*rhoinfd)/rhoinf)
    uinf2 = gammainf*(machcoef*machcoef*temp0)
! set the turbulent variables if transport variables are to be
! solved. we should be checking for rans equations here,
! however, this code is included in block res. the issue is
! that for frozen turbulence (or ank jacobian) we call the
! block_res with equationtype set to laminar even though we are
! actually solving the rans equations. the issue is that, the
! freestream turb variables will be changed to zero, thus
! changing the solution. insteady we check if nw > nwf which
! will accomplish the same thing.
    if (nw .gt. nwf) then
      nuinfd = (muinfd-muinf*rhoinfd/rhoinf)/rhoinf
      nuinf = muinf/rhoinf
      select case  (turbmodel) 
      case (spalartallmaras, spalartallmarasedwards) 
        winfd(itu1) = sanuknowneddyratio_d(eddyvisinfratio, nuinf, &
&         nuinfd, winf(itu1))
!=============================================================
      case (komegawilcox, komegamodified, mentersst) 
        winfd(itu1) = turbintensityinf**2*1.5_realtype*uinf2d
        winf(itu1) = 1.5_realtype*uinf2*turbintensityinf**2
        temp0 = winf(itu1)/(eddyvisinfratio*nuinf)
        winfd(itu2) = (winfd(itu1)-temp0*eddyvisinfratio*nuinfd)/(&
&         eddyvisinfratio*nuinf)
        winf(itu2) = temp0
!both are consistent with https://www.cfd-online.com/wiki/turbulence_free-stream_boundary_conditions,
! nasa https://turbmodels.larc.nasa.gov/sst.html has slightly different values
!the nasa ref specify that the freestream turbulent viscosity should be between 10-5 and 10-2 times freestream laminar viscosity.
! not clear why eddyvisinfratio default to 0.009
!this ref suggests similar things: k determined so that nutinf = nuinf * 0.009
! https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.901.7078&rep=rep1&type=pdf
!=============================================================
      case (ktau) 
        winfd(itu1) = turbintensityinf**2*1.5_realtype*uinf2d
        winf(itu1) = 1.5_realtype*uinf2*turbintensityinf**2
        temp0 = nuinf/winf(itu1)
        winfd(itu2) = eddyvisinfratio*(nuinfd-temp0*winfd(itu1))/winf(&
&         itu1)
        winf(itu2) = eddyvisinfratio*temp0
!=============================================================
      case (v2f) 
        winfd(itu1) = turbintensityinf**2*1.5_realtype*uinf2d
        winf(itu1) = 1.5_realtype*uinf2*turbintensityinf**2
        temp0 = winf(itu1)*winf(itu1)/(eddyvisinfratio*nuinf)
        winfd(itu2) = 0.09_realtype*(2*winf(itu1)*winfd(itu1)-temp0*&
&         eddyvisinfratio*nuinfd)/(eddyvisinfratio*nuinf)
        winf(itu2) = 0.09_realtype*temp0
        winfd(itu3) = 0.666666_realtype*winfd(itu1)
        winf(itu3) = 0.666666_realtype*winf(itu1)
        winfd(itu4) = 0.0_8
        winf(itu4) = 0.0_realtype
      end select
    end if
! set the value of pinfcorr. in case a k-equation is present
! add 2/3 times rho*k.
    pinfcorrd = pinfd
    pinfcorr = pinf
    if (kpresent) then
      pinfcorrd = pinfd + two*third*(winf(itu1)*rhoinfd+rhoinf*winfd(&
&       itu1))
      pinfcorr = pinf + two*third*rhoinf*winf(itu1)
    end if
! compute the free stream total energy.
    ktmp = zero
    if (kpresent) then
      ktmpd = winfd(itu1)
      ktmp = winf(itu1)
    else
      ktmpd = 0.0_8
    end if
    vinf = zero
    zinf = zero
    zinfd = 0.0_8
    vinfd = 0.0_8
    call etot_d(rhoinf, rhoinfd, uinf, uinfd, vinf, vinfd, zinf, zinfd, &
&         pinfcorr, pinfcorrd, ktmp, ktmpd, winf(irhoe), winfd(irhoe), &
&         kpresent)
  end subroutine referencestate_d

  subroutine referencestate()
!
!       the original version has been nuked since the computations are
!       no longer necessary when calling from python
!       this is the most compliclated routine in all of adflow. it is
!       stupidly complicated. this is most likely the reason your
!       derivatives are wrong. you don't understand this routine
!       and its effects.
!       this routine *requries* the following as input:
!       mach, pinfdim, tinfdim, rhoinfdim, rgasdim (machcoef non-sa
!        turbulence only)
!       optionally, pref, rhoref and tref are used if they are
!       are non-negative. this only happens when you want the equations
!       normalized by values other than the freestream
!      * this routine computes as output:
!      *   muinfdim, (unused anywhere in code)
!         pref, rhoref, tref, muref, timeref ('dimensional' reference)
!         pinf, pinfcorr, rhoinf, uinf, rgas, muinf, gammainf and winf
!         (non-dimensionalized values used in actual computations)
!
    use constants
    use paramturb
    use inputphysics, only : equations, mach, machcoef, musuthdim, &
&   tsuthdim, veldirfreestream, rgasdim, ssuthdim, eddyvisinfratio, &
&   turbmodel, turbintensityinf
    use flowvarrefstate, only : pinfdim, tinfdim, rhoinfdim, muinfdim,&
&   pref, rhoref, tref, muref, timeref, uref, href, pinf, pinfcorr, &
&   rhoinf, uinf, rgas, muinf, gammainf, winf, nw, nwf, kpresent, winf
    use flowutils_d, only : computegamma, etot
    use turbutils_d, only : sanuknowneddyratio
    implicit none
    integer(kind=inttype) :: sps, nn, mm, ierr
    real(kind=realtype) :: gm1, ratio
    real(kind=realtype) :: nuinf, ktmp, uinf2
    real(kind=realtype) :: vinf, zinf, tmp1(1), tmp2(1)
    intrinsic sqrt
    real(kind=realtype) :: arg1
    real(kind=realtype) :: result1
! compute the dimensional viscosity from sutherland's law
    muinfdim = musuthdim*((tsuthdim+ssuthdim)/(tinfdim+ssuthdim))*(&
&     tinfdim/tsuthdim)**1.5_realtype
! set the reference values. they *could* be different from the
! free-stream values for an internal flow simulation. for now,
! we just use the actual free stream values.
    pref = pinfdim
    tref = tinfdim
    rhoref = rhoinfdim
! compute the value of muref, such that the nondimensional
! equations are identical to the dimensional ones.
! note that in the non-dimensionalization of muref there is
! a reference length. however this reference length is 1.0
! in this code, because the coordinates are converted to
! meters.
    muref = sqrt(pref*rhoref)
! compute timeref for a correct nondimensionalization of the
! unsteady equations. some story as for the reference viscosity
! concerning the reference length.
    timeref = sqrt(rhoref/pref)
    href = pref/rhoref
    uref = sqrt(href)
! compute the nondimensional pressure, density, velocity,
! viscosity and gas constant.
    pinf = pinfdim/pref
    rhoinf = rhoinfdim/rhoref
    arg1 = gammainf*pinf/rhoinf
    result1 = sqrt(arg1)
    uinf = mach*result1
    rgas = rgasdim*rhoref*tref/pref
    muinf = muinfdim/muref
    tmp1(1) = tinfdim
    call computegamma(tmp1, tmp2, 1)
    gammainf = tmp2(1)
! ----------------------------------------
!      compute the final winf
! ----------------------------------------
! allocate the memory for winf if necessary
! zero out the winf first
    winf(:) = zero
! set the reference value of the flow variables, except the total
! energy. this will be computed at the end of this routine.
    winf(irho) = rhoinf
    winf(ivx) = uinf*veldirfreestream(1)
    winf(ivy) = uinf*veldirfreestream(2)
    winf(ivz) = uinf*veldirfreestream(3)
! compute the velocity squared based on machcoef. this gives a
! better indication of the 'speed' of the flow so the turubulence
! intensity ration is more meaningful especially for moving
! geometries. (not used in sa model)
    uinf2 = machcoef*machcoef*gammainf*pinf/rhoinf
! set the turbulent variables if transport variables are to be
! solved. we should be checking for rans equations here,
! however, this code is included in block res. the issue is
! that for frozen turbulence (or ank jacobian) we call the
! block_res with equationtype set to laminar even though we are
! actually solving the rans equations. the issue is that, the
! freestream turb variables will be changed to zero, thus
! changing the solution. insteady we check if nw > nwf which
! will accomplish the same thing.
    if (nw .gt. nwf) then
      nuinf = muinf/rhoinf
      select case  (turbmodel) 
      case (spalartallmaras, spalartallmarasedwards) 
        winf(itu1) = sanuknowneddyratio(eddyvisinfratio, nuinf)
!=============================================================
      case (komegawilcox, komegamodified, mentersst) 
        winf(itu1) = 1.5_realtype*uinf2*turbintensityinf**2
        winf(itu2) = winf(itu1)/(eddyvisinfratio*nuinf)
!both are consistent with https://www.cfd-online.com/wiki/turbulence_free-stream_boundary_conditions,
! nasa https://turbmodels.larc.nasa.gov/sst.html has slightly different values
!the nasa ref specify that the freestream turbulent viscosity should be between 10-5 and 10-2 times freestream laminar viscosity.
! not clear why eddyvisinfratio default to 0.009
!this ref suggests similar things: k determined so that nutinf = nuinf * 0.009
! https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.901.7078&rep=rep1&type=pdf
!=============================================================
      case (ktau) 
        winf(itu1) = 1.5_realtype*uinf2*turbintensityinf**2
        winf(itu2) = eddyvisinfratio*nuinf/winf(itu1)
!=============================================================
      case (v2f) 
        winf(itu1) = 1.5_realtype*uinf2*turbintensityinf**2
        winf(itu2) = 0.09_realtype*winf(itu1)**2/(eddyvisinfratio*nuinf)
        winf(itu3) = 0.666666_realtype*winf(itu1)
        winf(itu4) = 0.0_realtype
      end select
    end if
! set the value of pinfcorr. in case a k-equation is present
! add 2/3 times rho*k.
    pinfcorr = pinf
    if (kpresent) pinfcorr = pinf + two*third*rhoinf*winf(itu1)
! compute the free stream total energy.
    ktmp = zero
    if (kpresent) ktmp = winf(itu1)
    vinf = zero
    zinf = zero
    call etot(rhoinf, uinf, vinf, zinf, pinfcorr, ktmp, winf(irhoe), &
&       kpresent)
  end subroutine referencestate
! ----------------------------------------------------------------------
!                                                                      |
!                    no tapenade routine below this line               |
!                                                                      |
! ----------------------------------------------------------------------

end module initializeflow_d

