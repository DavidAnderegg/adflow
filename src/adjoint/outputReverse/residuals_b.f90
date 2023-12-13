!        generated by tapenade     (inria, ecuador team)
!  tapenade 3.16 (develop) - 10 nov 2023 18:24
!
module residuals_b
  implicit none

contains
  subroutine residual_block()
!
!       residual computes the residual of the mean flow equations on
!       the current mg level.
!
    use blockpointers
    use cgnsgrid
    use flowvarrefstate
    use inputiteration
    use inputdiscretization
    use inputtimespectral
! added by hdn
    use inputunsteady
    use iteration
    use inputadjoint
    use flowutils_b, only : computespeedofsoundsquared, &
&   allnodalgradients
    use fluxes_b
    implicit none
!
!      local variables.
!
    integer(kind=inttype) :: discr
    integer(kind=inttype) :: i, j, k, l
! for loops of ale
    integer(kind=inttype) :: iale, jale, kale, lale, male
    real(kind=realtype), parameter :: k1=1.05_realtype
! the line below is only used for the low-speed preconditioner part of this routine
! random given number
    real(kind=realtype), parameter :: k2=0.6_realtype
! mach number preconditioner activation
    real(kind=realtype), parameter :: m0=0.2_realtype
    real(kind=realtype), parameter :: alpha=0_realtype
    real(kind=realtype), parameter :: delta=0_realtype
!real(kind=realtype), parameter :: hinf = 2_realtype ! test phase
! test phase
    real(kind=realtype), parameter :: cpres=4.18_realtype
    real(kind=realtype), parameter :: temp=297.15_realtype
!
!     local variables
!
    real(kind=realtype) :: k3, h, velxrho, velyrho, velzrho, sos, hinf
    real(kind=realtype) :: resm, a11, a12, a13, a14, a15, a21, a22, a23&
&   , a24, a25, a31, a32, a33, a34, a35
    real(kind=realtype) :: a41, a42, a43, a44, a45, a51, a52, a53, a54, &
&   a55, b11, b12, b13, b14, b15
    real(kind=realtype) :: b21, b22, b23, b24, b25, b31, b32, b33, b34, &
&   b35
    real(kind=realtype) :: b41, b42, b43, b44, b45, b51, b52, b53, b54, &
&   b55
    real(kind=realtype) :: rhohdash, betamr2
    real(kind=realtype) :: g, q
    real(kind=realtype) :: b1, b2, b3, b4, b5
    real(kind=realtype) :: dwo(nwf)
    logical :: finegrid
    intrinsic abs
    intrinsic sqrt
    intrinsic max
    intrinsic min
    intrinsic real
    real(kind=realtype) :: x1
    real(realtype) :: x2
    real(realtype) :: x3
    real(kind=realtype) :: abs0
    real(realtype) :: max1
    real(realtype) :: max2
! set the value of rfil, which controls the fraction of the old
! dissipation residual to be used. this is only for the runge-kutta
! schemes; for other smoothers rfil is simply set to 1.0.
! note the index rkstage+1 for cdisrk. the reason is that the
! residual computation is performed before rkstage is incremented.
    if (smoother .eq. rungekutta) then
      rfil = cdisrk(rkstage+1)
    else
      rfil = one
    end if
! set the value of the discretization, depending on the grid level,
! and the logical finegrid, which indicates whether or not this
! is the finest grid level of the current mg cycle.
    discr = spacediscrcoarse
    if (currentlevel .eq. 1) discr = spacediscr
    finegrid = .false.
    if (currentlevel .eq. groundlevel) finegrid = .true.
! ===========================================================
!
! assuming ale has nothing to do with mg
! the geometric data will be interpolated if in md mode
!
! ===========================================================
! ===========================================================
!
! the fluxes are calculated as usual
!
! ===========================================================
    call inviscidcentralflux()
    select case  (discr) 
    case (dissscalar) 
! standard scalar dissipation scheme.
      if (finegrid) then
        if (.not.lumpeddiss) then
          call invisciddissfluxscalar()
        else
          call invisciddissfluxscalarapprox()
        end if
      end if
    case (dissmatrix) 
!===========================================================
! matrix dissipation scheme.
      if (finegrid) then
        if (.not.lumpeddiss) then
          call invisciddissfluxmatrix()
        else
          call invisciddissfluxmatrixapprox()
        end if
      end if
    case (upwind) 
!===========================================================
! dissipation via an upwind scheme.
      call inviscidupwindflux(finegrid)
    end select
!-------------------------------------------------------
! lastly, recover the old s[i,j,k], sface[i,j,k]
! this shall be done before difussive and source terms
! are computed.
!-------------------------------------------------------
    if (viscous) then
      if (rfil .ge. 0.) then
        abs0 = rfil
      else
        abs0 = -rfil
      end if
! only compute viscous fluxes if rfil > 0
      if (abs0 .gt. thresholdreal) then
! not lumpeddiss means it isn't the pc...call the vicousflux
        if (.not.lumpeddiss) then
          call computespeedofsoundsquared()
          call allnodalgradients()
          call viscousflux()
        else
! this is a pc calc...only include viscous fluxes if viscpc
! is used
! if full visc is true, also need full viscous terms, even if
! lumpeddiss is true
          call computespeedofsoundsquared()
          if (viscpc) then
            call allnodalgradients()
            call viscousflux()
          else
            call viscousfluxapprox()
          end if
        end if
      end if
    end if
!===========================================================
! add the dissipative and possibly viscous fluxes to the
! euler fluxes. loop over the owned cells and add fw to dw.
! also multiply by iblank so that no updates occur in holes
    if (lowspeedpreconditioner) then
      do k=2,kl
        do j=2,jl
          do i=2,il
!    compute speed of sound
            sos = sqrt(gamma(i, j, k)*p(i, j, k)/w(i, j, k, irho))
! compute velocities without rho from state vector
!      (w is pointer.. see type blocktype setup in block.f90)
!      w(0:ib,0:jb,0:kb,1:nw) is allocated in block.f90
!      these are per definition nw=[rho,u,v,w,rhoee]
!      so the velocity is simply just taken out below...
!      we do not have to divide with rho since it is already
!      without rho...
! ivx: l. 60 in constants.f90
            velxrho = w(i, j, k, ivx)
            velyrho = w(i, j, k, ivy)
            velzrho = w(i, j, k, ivz)
            q = velxrho**2 + velyrho**2 + velzrho**2
            resm = sqrt(q)/sos
! resm above is used as m_a (thesis) and m (paper 2015)
! and is the free stream mach number
! see routine setup above:
! l. 30: real(kind=realtype), parameter :: k1 = 1.05_realtype
! random given number for k2:
! l. 31: real(kind=realtype), parameter :: k2 = 0.6_realtype
! mach number preconditioner activation for k3:
! l. 32: real(kind=realtype), parameter :: m0 = 0.2_realtype
!
!    compute k3
! eq. 2.7 in garg 2015. k1, m0 and resm are scalars
!
! unfortunately, garg has switched the k1 and k3 here in the
! code. in both paper and thesis it is k3 that is used to det-
! ermine k1 below
!
!    compute k3
            k3 = k1*(1+(1-k1*m0**2)*resm**2/(k1*m0**4))
            if (k3*(velxrho**2+velyrho**2+velzrho**2) .lt. k2*(winf(ivx)&
&               **2+winf(ivy)**2+winf(ivz)**2)) then
              x1 = k2*(winf(ivx)**2+winf(ivy)**2+winf(ivz)**2)
            else
              x1 = k3*(velxrho**2+velyrho**2+velzrho**2)
            end if
            if (x1 .gt. sos**2) then
              betamr2 = sos**2
            else
              betamr2 = x1
            end if
! above, the winf is the free stream velocity
!
! should this first line's first element have sos^4 or sos^2
            a11 = betamr2*(1/sos**4)
            a12 = zero
            a13 = zero
            a14 = zero
            a15 = (-betamr2)/sos**4
            a21 = one*velxrho/sos**2
            a22 = one*w(i, j, k, irho)
            a23 = zero
            a24 = zero
            a25 = one*(-velxrho)/sos**2
            a31 = one*velyrho/sos**2
            a32 = zero
            a33 = one*w(i, j, k, irho)
            a34 = zero
            a35 = one*(-velyrho)/sos**2
            a41 = one*velzrho/sos**2
            a42 = zero
            a43 = zero
            a44 = one*w(i, j, k, irho)
            a45 = zero + one*(-velzrho)/sos**2
! mham: seems he fixed the above line an irregular way?
            a51 = one*(1/(gamma(i, j, k)-1)+resm**2/2)
            a52 = one*w(i, j, k, irho)*velxrho
            a53 = one*w(i, j, k, irho)*velyrho
            a54 = one*w(i, j, k, irho)*velzrho
            a55 = one*((-(resm**2))/2)
            b11 = a11*(gamma(i, j, k)-1)*q/2 + a12*(-velxrho)/w(i, j, k&
&             , irho) + a13*(-velyrho)/w(i, j, k, irho) + a14*(-velzrho)&
&             /w(i, j, k, irho) + a15*((gamma(i, j, k)-1)*q/2-sos**2)
            b12 = a11*(1-gamma(i, j, k))*velxrho + a12*1/w(i, j, k, irho&
&             ) + a15*(1-gamma(i, j, k))*velxrho
            b13 = a11*(1-gamma(i, j, k))*velyrho + a13/w(i, j, k, irho) &
&             + a15*(1-gamma(i, j, k))*velyrho
            b14 = a11*(1-gamma(i, j, k))*velzrho + a14/w(i, j, k, irho) &
&             + a15*(1-gamma(i, j, k))*velzrho
            b15 = a11*(gamma(i, j, k)-1) + a15*(gamma(i, j, k)-1)
            b21 = a21*(gamma(i, j, k)-1)*q/2 + a22*(-velxrho)/w(i, j, k&
&             , irho) + a23*(-velyrho)/w(i, j, k, irho) + a24*(-velzrho)&
&             /w(i, j, k, irho) + a25*((gamma(i, j, k)-1)*q/2-sos**2)
            b22 = a21*(1-gamma(i, j, k))*velxrho + a22/w(i, j, k, irho) &
&             + a25*(1-gamma(i, j, k))*velxrho
            b23 = a21*(1-gamma(i, j, k))*velyrho + a23*1/w(i, j, k, irho&
&             ) + a25*(1-gamma(i, j, k))*velyrho
            b24 = a21*(1-gamma(i, j, k))*velzrho + a24*1/w(i, j, k, irho&
&             ) + a25*(1-gamma(i, j, k))*velzrho
            b25 = a21*(gamma(i, j, k)-1) + a25*(gamma(i, j, k)-1)
            b31 = a31*(gamma(i, j, k)-1)*q/2 + a32*(-velxrho)/w(i, j, k&
&             , irho) + a33*(-velyrho)/w(i, j, k, irho) + a34*(-velzrho)&
&             /w(i, j, k, irho) + a35*((gamma(i, j, k)-1)*q/2-sos**2)
            b32 = a31*(1-gamma(i, j, k))*velxrho + a32/w(i, j, k, irho) &
&             + a35*(1-gamma(i, j, k))*velxrho
            b33 = a31*(1-gamma(i, j, k))*velyrho + a33*1/w(i, j, k, irho&
&             ) + a35*(1-gamma(i, j, k))*velyrho
            b34 = a31*(1-gamma(i, j, k))*velzrho + a34*1/w(i, j, k, irho&
&             ) + a35*(1-gamma(i, j, k))*velzrho
            b35 = a31*(gamma(i, j, k)-1) + a35*(gamma(i, j, k)-1)
            b41 = a41*(gamma(i, j, k)-1)*q/2 + a42*(-velxrho)/w(i, j, k&
&             , irho) + a43*(-velyrho)/w(i, j, k, irho) + a44*(-velzrho)&
&             /w(i, j, k, irho) + a45*((gamma(i, j, k)-1)*q/2-sos**2)
            b42 = a41*(1-gamma(i, j, k))*velxrho + a42/w(i, j, k, irho) &
&             + a45*(1-gamma(i, j, k))*velxrho
            b43 = a41*(1-gamma(i, j, k))*velyrho + a43*1/w(i, j, k, irho&
&             ) + a45*(1-gamma(i, j, k))*velyrho
            b44 = a41*(1-gamma(i, j, k))*velzrho + a44*1/w(i, j, k, irho&
&             ) + a45*(1-gamma(i, j, k))*velzrho
            b45 = a41*(gamma(i, j, k)-1) + a45*(gamma(i, j, k)-1)
            b51 = a51*(gamma(i, j, k)-1)*q/2 + a52*(-velxrho)/w(i, j, k&
&             , irho) + a53*(-velyrho)/w(i, j, k, irho) + a54*(-velzrho)&
&             /w(i, j, k, irho) + a55*((gamma(i, j, k)-1)*q/2-sos**2)
            b52 = a51*(1-gamma(i, j, k))*velxrho + a52/w(i, j, k, irho) &
&             + a55*(1-gamma(i, j, k))*velxrho
            b53 = a51*(1-gamma(i, j, k))*velyrho + a53*1/w(i, j, k, irho&
&             ) + a55*(1-gamma(i, j, k))*velyrho
            b54 = a51*(1-gamma(i, j, k))*velzrho + a54*1/w(i, j, k, irho&
&             ) + a55*(1-gamma(i, j, k))*velzrho
            b55 = a51*(gamma(i, j, k)-1) + a55*(gamma(i, j, k)-1)
! dwo is the orginal redisual
            do l=1,nwf
              x2 = real(iblank(i, j, k), realtype)
              if (x2 .lt. zero) then
                max1 = zero
              else
                max1 = x2
              end if
              dwo(l) = (dw(i, j, k, l)+fw(i, j, k, l))*max1
            end do
            dw(i, j, k, 1) = b11*dwo(1) + b12*dwo(2) + b13*dwo(3) + b14*&
&             dwo(4) + b15*dwo(5)
            dw(i, j, k, 2) = b21*dwo(1) + b22*dwo(2) + b23*dwo(3) + b24*&
&             dwo(4) + b25*dwo(5)
            dw(i, j, k, 3) = b31*dwo(1) + b32*dwo(2) + b33*dwo(3) + b34*&
&             dwo(4) + b35*dwo(5)
            dw(i, j, k, 4) = b41*dwo(1) + b42*dwo(2) + b43*dwo(3) + b44*&
&             dwo(4) + b45*dwo(5)
            dw(i, j, k, 5) = b51*dwo(1) + b52*dwo(2) + b53*dwo(3) + b54*&
&             dwo(4) + b55*dwo(5)
          end do
        end do
! end of lowspeedpreconditioners three cells loops

      end do
    else
! else.. i.e. if we do not have preconditioner turned on...
      do l=1,nwf
        do k=2,kl
          do j=2,jl
            do i=2,il
              x3 = real(iblank(i, j, k), realtype)
              if (x3 .lt. zero) then
                max2 = zero
              else
                max2 = x3
              end if
              dw(i, j, k, l) = (dw(i, j, k, l)+fw(i, j, k, l))*max2
            end do
          end do
        end do
      end do
    end if
  end subroutine residual_block

!  differentiation of sourceterms_block in reverse (adjoint) mode (with options noisize i4 dr8 r8):
!   gradient     of useful results: uref pref *w *dw *vol actuatorregions.force
!                actuatorregions.heat actuatorregions.volume plocal
!   with respect to varying inputs: uref pref *w *dw *vol actuatorregions.force
!                actuatorregions.heat actuatorregions.volume plocal
!   rw status of diff variables: uref:incr pref:incr *w:incr *dw:in-out
!                *vol:incr actuatorregions.force:incr actuatorregions.heat:incr
!                actuatorregions.volume:incr plocal:in-out
!   plus diff mem management of: w:in dw:in vol:in
  subroutine sourceterms_block_b0(nn, res, iregion, plocal, plocald)
! apply the source terms for the given block. assume that the
! block pointers are already set.
    use constants
    use actuatorregiondata
    use blockpointers, only : vol, vold, dw, dwd, w, wd
    use flowvarrefstate, only : pref, prefd, uref, urefd, lref
    use communication
    use iteration, only : ordersconverged
    implicit none
! input
    integer(kind=inttype), intent(in) :: nn, iregion
    logical, intent(in) :: res
    real(kind=realtype), intent(inout) :: plocal
    real(kind=realtype), intent(inout) :: plocald
! working
    integer(kind=inttype) :: i, j, k, ii, istart, iend
    real(kind=realtype) :: ftmp(3), vx, vy, vz, f_fact(3), q_fact, qtmp&
&   , redim, factor, ostart, oend
    real(kind=realtype) :: ftmpd(3), vxd, vyd, vzd, f_factd(3), q_factd&
&   , qtmpd, redimd
    real(kind=realtype) :: thrustt(3), thrusttnorm
    real(kind=realtype), dimension(3) :: tempd
    real(kind=realtype) :: temp
    real(kind=realtype) :: tempd0
    real(kind=realtype) :: temp0
    real(kind=realtype) :: tempd1
    real(kind=realtype) :: tempd2
    integer :: branch
    redim = pref*uref
! compute the relaxation factor based on the ordersconverged
! how far we are into the ramp:
    if (ordersconverged .lt. actuatorregions(iregion)%relaxstart) then
      call pushcontrol1b(0)
      factor = zero
    else if (ordersconverged .gt. actuatorregions(iregion)%relaxend) &
&   then
      call pushcontrol1b(1)
      factor = one
    else
      call pushcontrol1b(1)
! in between
      ostart = actuatorregions(iregion)%relaxstart
      oend = actuatorregions(iregion)%relaxend
      factor = (ordersconverged-ostart)/(oend-ostart)
    end if
! if using the uniform force distribution
    if (actuatorregions(iregion)%acttype .eq. 'uniform') then
! compute the constant force factor
      f_fact = factor*actuatorregions(iregion)%force/actuatorregions(&
&       iregion)%volume/pref
      call pushcontrol1b(0)
    else
      call pushcontrol1b(1)
    end if
! heat factor. this is heat added per unit volume per unit time
    q_fact = factor*actuatorregions(iregion)%heat/actuatorregions(&
&     iregion)%volume/(pref*uref*lref*lref)
! loop over the ranges for this block
    istart = actuatorregions(iregion)%blkptr(nn-1) + 1
    iend = actuatorregions(iregion)%blkptr(nn)
! if using the uniform force distribution
    if (actuatorregions(iregion)%acttype .eq. 'uniform') then
      call pushcontrol1b(0)
    else
      call pushcontrol1b(1)
    end if
! if using the simple propeller force distribution
    if (actuatorregions(iregion)%acttype .eq. 'simpleprop') then
      q_factd = 0.0_8
      redimd = 0.0_8
!$bwd-of ii-loop 
      do ii=istart,iend
! extract the cell id.
        i = actuatorregions(iregion)%cellids(1, ii)
        j = actuatorregions(iregion)%cellids(2, ii)
        k = actuatorregions(iregion)%cellids(3, ii)
!thrustt = actuatorregions(iregion)%force
!thrusttnorm = sqrt((thrustt(1)**2 + thrustt(2)**2 + thrustt(3)**2))
        ftmp = factor*actuatorregions(iregion)%thrustvec(:, ii)*vol(i, j&
&         , k)*actuatorregions(iregion)%thrust/pref
        ftmp = ftmp + factor*actuatorregions(iregion)%swirlvec(:, ii)*&
&         actuatorregions(iregion)%thrust/pref
        vx = w(i, j, k, ivx)
        vy = w(i, j, k, ivy)
        vz = w(i, j, k, ivz)
        if (res) then
          ftmpd = 0.0_8
          ftmpd(1) = ftmpd(1) - vx*dwd(i, j, k, irhoe)
          vxd = -(ftmp(1)*dwd(i, j, k, irhoe))
          ftmpd(2) = ftmpd(2) - vy*dwd(i, j, k, irhoe)
          vyd = -(ftmp(2)*dwd(i, j, k, irhoe))
          ftmpd(3) = ftmpd(3) - vz*dwd(i, j, k, irhoe)
          vzd = -(ftmp(3)*dwd(i, j, k, irhoe))
          qtmpd = -dwd(i, j, k, irhoe)
          ftmpd = ftmpd - dwd(i, j, k, imx:imz)
        else
          ftmpd = 0.0_8
          tempd1 = redim*plocald
          redimd = redimd + (vx*ftmp(1)+vy*ftmp(2)+vz*ftmp(3))*plocald
          vxd = ftmp(1)*tempd1
          ftmpd(1) = ftmpd(1) + vx*tempd1
          vyd = ftmp(2)*tempd1
          ftmpd(2) = ftmpd(2) + vy*tempd1
          vzd = ftmp(3)*tempd1
          ftmpd(3) = ftmpd(3) + vz*tempd1
          qtmpd = 0.0_8
        end if
        tempd1 = sum(actuatorregions(iregion)%thrustvec(:, ii)*ftmpd)*&
&         factor*actuatorregions(iregion)%thrust/pref
        wd(i, j, k, ivz) = wd(i, j, k, ivz) + vzd
        wd(i, j, k, ivy) = wd(i, j, k, ivy) + vyd
        wd(i, j, k, ivx) = wd(i, j, k, ivx) + vxd
        vold(i, j, k) = vold(i, j, k) + q_fact*qtmpd + tempd1
        q_factd = q_factd + vol(i, j, k)*qtmpd
        prefd = prefd - factor*actuatorregions(iregion)%thrust*sum(&
&         actuatorregions(iregion)%swirlvec(:, ii)*ftmpd)/pref**2 - vol(&
&         i, j, k)*tempd1/pref
      end do
    else
      q_factd = 0.0_8
      redimd = 0.0_8
    end if
    call popcontrol1b(branch)
    if (branch .eq. 0) then
      f_factd = 0.0_8
!$bwd-of ii-loop 
      do ii=istart,iend
! extract the cell id.
        i = actuatorregions(iregion)%cellids(1, ii)
        j = actuatorregions(iregion)%cellids(2, ii)
        k = actuatorregions(iregion)%cellids(3, ii)
! this actually gets the force
        ftmp = vol(i, j, k)*f_fact
        vx = w(i, j, k, ivx)
        vy = w(i, j, k, ivy)
        vz = w(i, j, k, ivz)
! this gets the heat addition rate
        if (res) then
          ftmpd = 0.0_8
          ftmpd(1) = ftmpd(1) - vx*dwd(i, j, k, irhoe)
          vxd = -(ftmp(1)*dwd(i, j, k, irhoe))
          ftmpd(2) = ftmpd(2) - vy*dwd(i, j, k, irhoe)
          vyd = -(ftmp(2)*dwd(i, j, k, irhoe))
          ftmpd(3) = ftmpd(3) - vz*dwd(i, j, k, irhoe)
          vzd = -(ftmp(3)*dwd(i, j, k, irhoe))
          qtmpd = -dwd(i, j, k, irhoe)
          ftmpd = ftmpd - dwd(i, j, k, imx:imz)
        else
          ftmpd = 0.0_8
          tempd1 = redim*plocald
          redimd = redimd + (vx*ftmp(1)+vy*ftmp(2)+vz*ftmp(3))*plocald
          vxd = ftmp(1)*tempd1
          ftmpd(1) = ftmpd(1) + vx*tempd1
          vyd = ftmp(2)*tempd1
          ftmpd(2) = ftmpd(2) + vy*tempd1
          vzd = ftmp(3)*tempd1
          ftmpd(3) = ftmpd(3) + vz*tempd1
          qtmpd = 0.0_8
        end if
        vold(i, j, k) = vold(i, j, k) + q_fact*qtmpd + sum(f_fact*ftmpd)
        q_factd = q_factd + vol(i, j, k)*qtmpd
        wd(i, j, k, ivz) = wd(i, j, k, ivz) + vzd
        wd(i, j, k, ivy) = wd(i, j, k, ivy) + vyd
        wd(i, j, k, ivx) = wd(i, j, k, ivx) + vxd
        f_factd = f_factd + vol(i, j, k)*ftmpd
      end do
    else
      f_factd = 0.0_8
    end if
    temp = lref*lref*actuatorregions(iregion)%volume
    temp0 = temp*pref*uref
    tempd1 = factor*q_factd/temp0
    actuatorregionsd(iregion)%heat = actuatorregionsd(iregion)%heat + &
&     tempd1
    tempd2 = -(actuatorregions(iregion)%heat*tempd1/temp0)
    actuatorregionsd(iregion)%volume = actuatorregionsd(iregion)%volume &
&     + lref**2*pref*uref*tempd2
    prefd = prefd + uref*temp*tempd2
    urefd = urefd + pref*temp*tempd2
    call popcontrol1b(branch)
    if (branch .eq. 0) then
      tempd = factor*f_factd/(actuatorregions(iregion)%volume*pref)
      actuatorregionsd(iregion)%force = actuatorregionsd(iregion)%force &
&       + tempd
      tempd0 = -(sum(actuatorregions(iregion)%force*tempd)/(&
&       actuatorregions(iregion)%volume*pref))
      actuatorregionsd(iregion)%volume = actuatorregionsd(iregion)%&
&       volume + pref*tempd0
      prefd = prefd + actuatorregions(iregion)%volume*tempd0
    end if
    call popcontrol1b(branch)
    prefd = prefd + uref*redimd
    urefd = urefd + pref*redimd
  end subroutine sourceterms_block_b0

!  differentiation of sourceterms_block in reverse (adjoint) mode (with options noisize i4 dr8 r8):
!   gradient     of useful results: uref pref *w *dw actuatorregions.force
!                actuatorregions.thrust actuatorregions.heat plocal
!   with respect to varying inputs: uref pref *w *dw actuatorregions.force
!                actuatorregions.thrust actuatorregions.heat plocal
!   rw status of diff variables: uref:incr pref:incr *w:incr *dw:in-out
!                actuatorregions.force:incr actuatorregions.thrust:incr
!                actuatorregions.heat:incr plocal:in-out
!   plus diff mem management of: w:in dw:in
  subroutine sourceterms_block_b(nn, res, iregion, plocal, plocald)
! apply the source terms for the given block. assume that the
! block pointers are already set.
    use constants
    use actuatorregiondata
    use blockpointers, only : vol, vold, dw, dwd, w, wd
    use flowvarrefstate, only : pref, prefd, uref, urefd, lref
    use communication
    use iteration, only : ordersconverged
    implicit none
! input
    integer(kind=inttype), intent(in) :: nn, iregion
    logical, intent(in) :: res
    real(kind=realtype), intent(inout) :: plocal
    real(kind=realtype), intent(inout) :: plocald
! working
    integer(kind=inttype) :: i, j, k, ii, istart, iend
    real(kind=realtype) :: ftmp(3), vx, vy, vz, f_fact(3), q_fact, qtmp&
&   , redim, factor, ostart, oend
    real(kind=realtype) :: ftmpd(3), vxd, vyd, vzd, f_factd(3), q_factd&
&   , qtmpd, redimd
    real(kind=realtype) :: thrustt(3), thrusttnorm
    real(kind=realtype), dimension(3) :: tempd
    real(kind=realtype) :: temp
    real(kind=realtype) :: temp0
    real(kind=realtype) :: tempd0
    real(kind=realtype) :: tempd1
    integer :: branch
    redim = pref*uref
! compute the relaxation factor based on the ordersconverged
! how far we are into the ramp:
    if (ordersconverged .lt. actuatorregions(iregion)%relaxstart) then
      call pushcontrol1b(0)
      factor = zero
    else if (ordersconverged .gt. actuatorregions(iregion)%relaxend) &
&   then
      call pushcontrol1b(1)
      factor = one
    else
      call pushcontrol1b(1)
! in between
      ostart = actuatorregions(iregion)%relaxstart
      oend = actuatorregions(iregion)%relaxend
      factor = (ordersconverged-ostart)/(oend-ostart)
    end if
! if using the uniform force distribution
    if (actuatorregions(iregion)%acttype .eq. 'uniform') then
! compute the constant force factor
      f_fact = factor*actuatorregions(iregion)%force/actuatorregions(&
&       iregion)%volume/pref
      call pushcontrol1b(0)
    else
      call pushcontrol1b(1)
    end if
! heat factor. this is heat added per unit volume per unit time
! loop over the ranges for this block
    istart = actuatorregions(iregion)%blkptr(nn-1) + 1
    iend = actuatorregions(iregion)%blkptr(nn)
! if using the uniform force distribution
    if (actuatorregions(iregion)%acttype .eq. 'uniform') then
      call pushcontrol1b(0)
    else
      call pushcontrol1b(1)
    end if
! if using the simple propeller force distribution
    if (actuatorregions(iregion)%acttype .eq. 'simpleprop') then
      q_factd = 0.0_8
      redimd = 0.0_8
!$bwd-of ii-loop 
      do ii=istart,iend
! extract the cell id.
        i = actuatorregions(iregion)%cellids(1, ii)
        j = actuatorregions(iregion)%cellids(2, ii)
        k = actuatorregions(iregion)%cellids(3, ii)
!thrustt = actuatorregions(iregion)%force
!thrusttnorm = sqrt((thrustt(1)**2 + thrustt(2)**2 + thrustt(3)**2))
        ftmp = factor*actuatorregions(iregion)%thrustvec(:, ii)*vol(i, j&
&         , k)*actuatorregions(iregion)%thrust/pref
        ftmp = ftmp + factor*actuatorregions(iregion)%swirlvec(:, ii)*&
&         actuatorregions(iregion)%thrust/pref
        vx = w(i, j, k, ivx)
        vy = w(i, j, k, ivy)
        vz = w(i, j, k, ivz)
        if (res) then
          ftmpd = 0.0_8
          ftmpd(1) = ftmpd(1) - vx*dwd(i, j, k, irhoe)
          vxd = -(ftmp(1)*dwd(i, j, k, irhoe))
          ftmpd(2) = ftmpd(2) - vy*dwd(i, j, k, irhoe)
          vyd = -(ftmp(2)*dwd(i, j, k, irhoe))
          ftmpd(3) = ftmpd(3) - vz*dwd(i, j, k, irhoe)
          vzd = -(ftmp(3)*dwd(i, j, k, irhoe))
          qtmpd = -dwd(i, j, k, irhoe)
          ftmpd = ftmpd - dwd(i, j, k, imx:imz)
        else
          ftmpd = 0.0_8
          tempd1 = redim*plocald
          redimd = redimd + (vx*ftmp(1)+vy*ftmp(2)+vz*ftmp(3))*plocald
          vxd = ftmp(1)*tempd1
          ftmpd(1) = ftmpd(1) + vx*tempd1
          vyd = ftmp(2)*tempd1
          ftmpd(2) = ftmpd(2) + vy*tempd1
          vzd = ftmp(3)*tempd1
          ftmpd(3) = ftmpd(3) + vz*tempd1
          qtmpd = 0.0_8
        end if
        wd(i, j, k, ivz) = wd(i, j, k, ivz) + vzd
        wd(i, j, k, ivy) = wd(i, j, k, ivy) + vyd
        wd(i, j, k, ivx) = wd(i, j, k, ivx) + vxd
        q_factd = q_factd + vol(i, j, k)*qtmpd
        tempd1 = sum(actuatorregions(iregion)%swirlvec(:, ii)*ftmpd)*&
&         factor/pref
        actuatorregionsd(iregion)%thrust = actuatorregionsd(iregion)%&
&         thrust + tempd1
        prefd = prefd - actuatorregions(iregion)%thrust*tempd1/pref
        tempd1 = vol(i, j, k)*sum(actuatorregions(iregion)%thrustvec(:, &
&         ii)*ftmpd)*factor/pref
        actuatorregionsd(iregion)%thrust = actuatorregionsd(iregion)%&
&         thrust + tempd1
        prefd = prefd - actuatorregions(iregion)%thrust*tempd1/pref
      end do
    else
      q_factd = 0.0_8
      redimd = 0.0_8
    end if
    call popcontrol1b(branch)
    if (branch .eq. 0) then
      f_factd = 0.0_8
!$bwd-of ii-loop 
      do ii=istart,iend
! extract the cell id.
        i = actuatorregions(iregion)%cellids(1, ii)
        j = actuatorregions(iregion)%cellids(2, ii)
        k = actuatorregions(iregion)%cellids(3, ii)
! this actually gets the force
        ftmp = vol(i, j, k)*f_fact
        vx = w(i, j, k, ivx)
        vy = w(i, j, k, ivy)
        vz = w(i, j, k, ivz)
! this gets the heat addition rate
        if (res) then
          ftmpd = 0.0_8
          ftmpd(1) = ftmpd(1) - vx*dwd(i, j, k, irhoe)
          vxd = -(ftmp(1)*dwd(i, j, k, irhoe))
          ftmpd(2) = ftmpd(2) - vy*dwd(i, j, k, irhoe)
          vyd = -(ftmp(2)*dwd(i, j, k, irhoe))
          ftmpd(3) = ftmpd(3) - vz*dwd(i, j, k, irhoe)
          vzd = -(ftmp(3)*dwd(i, j, k, irhoe))
          qtmpd = -dwd(i, j, k, irhoe)
          ftmpd = ftmpd - dwd(i, j, k, imx:imz)
        else
          ftmpd = 0.0_8
          tempd1 = redim*plocald
          redimd = redimd + (vx*ftmp(1)+vy*ftmp(2)+vz*ftmp(3))*plocald
          vxd = ftmp(1)*tempd1
          ftmpd(1) = ftmpd(1) + vx*tempd1
          vyd = ftmp(2)*tempd1
          ftmpd(2) = ftmpd(2) + vy*tempd1
          vzd = ftmp(3)*tempd1
          ftmpd(3) = ftmpd(3) + vz*tempd1
          qtmpd = 0.0_8
        end if
        q_factd = q_factd + vol(i, j, k)*qtmpd
        wd(i, j, k, ivz) = wd(i, j, k, ivz) + vzd
        wd(i, j, k, ivy) = wd(i, j, k, ivy) + vyd
        wd(i, j, k, ivx) = wd(i, j, k, ivx) + vxd
        f_factd = f_factd + vol(i, j, k)*ftmpd
      end do
    else
      f_factd = 0.0_8
    end if
    temp = actuatorregions(iregion)%volume*(lref*lref)
    temp0 = temp*pref*uref
    tempd0 = factor*q_factd/temp0
    actuatorregionsd(iregion)%heat = actuatorregionsd(iregion)%heat + &
&     tempd0
    tempd1 = -(temp*actuatorregions(iregion)%heat*tempd0/temp0)
    prefd = prefd + uref*tempd1
    urefd = urefd + pref*tempd1
    call popcontrol1b(branch)
    if (branch .eq. 0) then
      tempd = factor*f_factd/(actuatorregions(iregion)%volume*pref)
      actuatorregionsd(iregion)%force = actuatorregionsd(iregion)%force &
&       + tempd
      prefd = prefd - sum(actuatorregions(iregion)%force*tempd)/pref
    end if
    call popcontrol1b(branch)
    prefd = prefd + uref*redimd
    urefd = urefd + pref*redimd
  end subroutine sourceterms_block_b

  subroutine sourceterms_block(nn, res, iregion, plocal)
! apply the source terms for the given block. assume that the
! block pointers are already set.
    use constants
    use actuatorregiondata
    use blockpointers, only : vol, dw, w
    use flowvarrefstate, only : pref, uref, lref
    use communication
    use iteration, only : ordersconverged
    implicit none
! input
    integer(kind=inttype), intent(in) :: nn, iregion
    logical, intent(in) :: res
    real(kind=realtype), intent(inout) :: plocal
! working
    integer(kind=inttype) :: i, j, k, ii, istart, iend
    real(kind=realtype) :: ftmp(3), vx, vy, vz, f_fact(3), q_fact, qtmp&
&   , redim, factor, ostart, oend
    real(kind=realtype) :: thrustt(3), thrusttnorm
    redim = pref*uref
! compute the relaxation factor based on the ordersconverged
! how far we are into the ramp:
    if (ordersconverged .lt. actuatorregions(iregion)%relaxstart) then
      factor = zero
    else if (ordersconverged .gt. actuatorregions(iregion)%relaxend) &
&   then
      factor = one
    else
! in between
      ostart = actuatorregions(iregion)%relaxstart
      oend = actuatorregions(iregion)%relaxend
      factor = (ordersconverged-ostart)/(oend-ostart)
    end if
! if using the uniform force distribution
    if (actuatorregions(iregion)%acttype .eq. 'uniform') then
! compute the constant force factor
      f_fact = factor*actuatorregions(iregion)%force/actuatorregions(&
&       iregion)%volume/pref
    end if
! heat factor. this is heat added per unit volume per unit time
    q_fact = factor*actuatorregions(iregion)%heat/actuatorregions(&
&     iregion)%volume/(pref*uref*lref*lref)
! loop over the ranges for this block
    istart = actuatorregions(iregion)%blkptr(nn-1) + 1
    iend = actuatorregions(iregion)%blkptr(nn)
! if using the uniform force distribution
    if (actuatorregions(iregion)%acttype .eq. 'uniform') then
!$ad ii-loop
      do ii=istart,iend
! extract the cell id.
        i = actuatorregions(iregion)%cellids(1, ii)
        j = actuatorregions(iregion)%cellids(2, ii)
        k = actuatorregions(iregion)%cellids(3, ii)
! this actually gets the force
        ftmp = vol(i, j, k)*f_fact
        vx = w(i, j, k, ivx)
        vy = w(i, j, k, ivy)
        vz = w(i, j, k, ivz)
! this gets the heat addition rate
        qtmp = vol(i, j, k)*q_fact
        if (res) then
! momentum residuals
          dw(i, j, k, imx:imz) = dw(i, j, k, imx:imz) - ftmp
! energy residuals
          dw(i, j, k, irhoe) = dw(i, j, k, irhoe) - ftmp(1)*vx - ftmp(2)&
&           *vy - ftmp(3)*vz - qtmp
        else
! add in the local power contribution:
          plocal = plocal + (vx*ftmp(1)+vy*ftmp(2)+vz*ftmp(3))*redim
        end if
      end do
    end if
! if using the simple propeller force distribution
    if (actuatorregions(iregion)%acttype .eq. 'simpleprop') then
!$ad ii-loop
      do ii=istart,iend
! extract the cell id.
        i = actuatorregions(iregion)%cellids(1, ii)
        j = actuatorregions(iregion)%cellids(2, ii)
        k = actuatorregions(iregion)%cellids(3, ii)
!thrustt = actuatorregions(iregion)%force
!thrusttnorm = sqrt((thrustt(1)**2 + thrustt(2)**2 + thrustt(3)**2))
        ftmp = factor*actuatorregions(iregion)%thrustvec(:, ii)*vol(i, j&
&         , k)*actuatorregions(iregion)%thrust/pref
        ftmp = ftmp + factor*actuatorregions(iregion)%swirlvec(:, ii)*&
&         actuatorregions(iregion)%thrust/pref
        qtmp = vol(i, j, k)*q_fact
        vx = w(i, j, k, ivx)
        vy = w(i, j, k, ivy)
        vz = w(i, j, k, ivz)
        if (res) then
! momentum residuals
          dw(i, j, k, imx:imz) = dw(i, j, k, imx:imz) - ftmp
! energy residuals
          dw(i, j, k, irhoe) = dw(i, j, k, irhoe) - ftmp(1)*vx - ftmp(2)&
&           *vy - ftmp(3)*vz - qtmp
        else
! add in the local power contribution:
          plocal = plocal + (vx*ftmp(1)+vy*ftmp(2)+vz*ftmp(3))*redim
        end if
      end do
    end if
  end subroutine sourceterms_block

!  differentiation of initres_block in reverse (adjoint) mode (with options noisize i4 dr8 r8):
!   gradient     of useful results: *dw
!   with respect to varying inputs: *dw
!   rw status of diff variables: *dw:in-out
!   plus diff mem management of: dw:in
  subroutine initres_block_b(varstart, varend, nn, sps)
!
!       initres initializes the given range of the residual. either to
!       zero, steady computation, or to an unsteady term for the time
!       spectral and unsteady modes. for the coarser grid levels the
!       residual forcing term is taken into account.
!
    use blockpointers
    use flowvarrefstate
    use inputiteration
    use inputphysics
    use inputtimespectral
    use inputunsteady
    use iteration
    implicit none
!
!      subroutine arguments.
!
    integer(kind=inttype), intent(in) :: varstart, varend, nn, sps
!
!      local variables.
!
    integer(kind=inttype) :: mm, ll, ii, jj, i, j, k, l, m
    real(kind=realtype) :: oneoverdt, tmp
    real(kind=realtype), dimension(:, :, :, :), pointer :: ww, wsp, wsp1
    real(kind=realtype), dimension(:, :, :), pointer :: volsp
    integer :: branch
! return immediately of no variables are in the range.
    if (varend .ge. varstart) then
! determine the equation mode and act accordingly.
      select case  (equationmode) 
      case (steady) 
! steady state computation.
! determine the currently active multigrid level.
        if (currentlevel .eq. groundlevel) then
          call pushcontrol2b(1)
        else
          call pushcontrol2b(0)
        end if
      case default
        call pushcontrol2b(2)
      end select
      do l=varend,varstart,-1
        do j=jl,2,-1
          do i=il,2,-1
            dwd(i, j, kb, l) = 0.0_8
            dwd(i, j, ke, l) = 0.0_8
            dwd(i, j, 1, l) = 0.0_8
            dwd(i, j, 0, l) = 0.0_8
          end do
        end do
        do k=kb,0,-1
          do i=il,2,-1
            dwd(i, jb, k, l) = 0.0_8
            dwd(i, je, k, l) = 0.0_8
            dwd(i, 1, k, l) = 0.0_8
            dwd(i, 0, k, l) = 0.0_8
          end do
        end do
        do k=kb,0,-1
          do j=jb,0,-1
            dwd(ib, j, k, l) = 0.0_8
            dwd(ie, j, k, l) = 0.0_8
            dwd(1, j, k, l) = 0.0_8
            dwd(0, j, k, l) = 0.0_8
          end do
        end do
      end do
      call popcontrol2b(branch)
      if (branch .eq. 0) then
        do l=varend,varstart,-1
          do k=kl,2,-1
            do j=jl,2,-1
              do i=il,2,-1
                dwd(i, j, k, l) = 0.0_8
              end do
            end do
          end do
        end do
      else if (branch .eq. 1) then
        do l=varend,varstart,-1
          do k=kl,2,-1
            do j=jl,2,-1
              do i=il,2,-1
                dwd(i, j, k, l) = 0.0_8
              end do
            end do
          end do
        end do
      end if
    end if
  end subroutine initres_block_b

  subroutine initres_block(varstart, varend, nn, sps)
!
!       initres initializes the given range of the residual. either to
!       zero, steady computation, or to an unsteady term for the time
!       spectral and unsteady modes. for the coarser grid levels the
!       residual forcing term is taken into account.
!
    use blockpointers
    use flowvarrefstate
    use inputiteration
    use inputphysics
    use inputtimespectral
    use inputunsteady
    use iteration
    implicit none
!
!      subroutine arguments.
!
    integer(kind=inttype), intent(in) :: varstart, varend, nn, sps
!
!      local variables.
!
    integer(kind=inttype) :: mm, ll, ii, jj, i, j, k, l, m
    real(kind=realtype) :: oneoverdt, tmp
    real(kind=realtype), dimension(:, :, :, :), pointer :: ww, wsp, wsp1
    real(kind=realtype), dimension(:, :, :), pointer :: volsp
! return immediately of no variables are in the range.
    if (varend .lt. varstart) then
      return
    else
! determine the equation mode and act accordingly.
      select case  (equationmode) 
      case (steady) 
! steady state computation.
! determine the currently active multigrid level.
        if (currentlevel .eq. groundlevel) then
! ground level of the multigrid cycle. initialize the
! owned residuals to zero.
          do l=varstart,varend
            do k=2,kl
              do j=2,jl
                do i=2,il
                  dw(i, j, k, l) = zero
                end do
              end do
            end do
          end do
        else
! coarse grid level. initialize the owned cells to the
! residual forcing terms.
          do l=varstart,varend
            do k=2,kl
              do j=2,jl
                do i=2,il
                  dw(i, j, k, l) = wr(i, j, k, l)
                end do
              end do
            end do
          end do
        end if
      end select
! set the residual in the halo cells to zero. this is just
! to avoid possible problems. their values do not matter.
      do l=varstart,varend
        do k=0,kb
          do j=0,jb
            dw(0, j, k, l) = zero
            dw(1, j, k, l) = zero
            dw(ie, j, k, l) = zero
            dw(ib, j, k, l) = zero
          end do
        end do
        do k=0,kb
          do i=2,il
            dw(i, 0, k, l) = zero
            dw(i, 1, k, l) = zero
            dw(i, je, k, l) = zero
            dw(i, jb, k, l) = zero
          end do
        end do
        do j=2,jl
          do i=2,il
            dw(i, j, 0, l) = zero
            dw(i, j, 1, l) = zero
            dw(i, j, ke, l) = zero
            dw(i, j, kb, l) = zero
          end do
        end do
      end do
    end if
  end subroutine initres_block
! ----------------------------------------------------------------------
!                                                                      |
!                    no tapenade routine below this line               |
!                                                                      |
! ----------------------------------------------------------------------

end module residuals_b

