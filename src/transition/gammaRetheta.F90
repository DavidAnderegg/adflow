module GammaRethetaModel

    use constants

contains

    subroutine solve_local_Re_thetat_eq(Re_thetat_eq, i, j, k)
        use blockPointers
        use constants
        use paramTurb
        implicit None

        ! input/output variables
        integer(kind=intType), intent(in) :: i, j, k
        real(kind=realType), intent(out) :: Re_thetat_eq

        ! local variables
        real(kind=realType) :: U, U_inv, fact, dU_dx, dU_dy, dU_dz, dU_ds, Tu, F1, F2, F3, F 
        real(kind=realType) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
        real(kind=realType) :: lambda, thetat, residum, thetat_old, residum_old, thetat_new
        real(kind=realType) :: Re_thetat_eq_1, Re_thetat_eq_2
        integer(kind=intType) :: n


        U = sqrt(w(i, j, k, ivx)**2 + w(i, j, k, ivy)**2 + w(i, j, k, ivz)**2)
        U_inv = 1.0/U

        ! Compute the gradient of u in the cell center. Use is made
        ! of the fact that the surrounding normals sum up to zero,
        ! such that the cell i,j,k does not give a contribution.
        ! Since the gradient is scaled by a factor of 2*vol, we need to account for that

        fact = 1.0 / (vol(i, j, k) * 2.0)

        dudx = (w(i + 1, j, k, ivx) * si(i, j, k, 1) - w(i - 1, j, k, ivx) * si(i - 1, j, k, 1) &
              + w(i, j + 1, k, ivx) * sj(i, j, k, 1) - w(i, j - 1, k, ivx) * sj(i, j - 1, k, 1) &
              + w(i, j, k + 1, ivx) * sk(i, j, k, 1) - w(i, j, k - 1, ivx) * sk(i, j, k - 1, 1)) * fact
        dudy = (w(i + 1, j, k, ivx) * si(i, j, k, 2) - w(i - 1, j, k, ivx) * si(i - 1, j, k, 2) &
              + w(i, j + 1, k, ivx) * sj(i, j, k, 2) - w(i, j - 1, k, ivx) * sj(i, j - 1, k, 2) &
              + w(i, j, k + 1, ivx) * sk(i, j, k, 2) - w(i, j, k - 1, ivx) * sk(i, j, k - 1, 2)) * fact
        dudz = (w(i + 1, j, k, ivx) * si(i, j, k, 3) - w(i - 1, j, k, ivx) * si(i - 1, j, k, 3) &
              + w(i, j + 1, k, ivx) * sj(i, j, k, 3) - w(i, j - 1, k, ivx) * sj(i, j - 1, k, 3) &
              + w(i, j, k + 1, ivx) * sk(i, j, k, 3) - w(i, j, k - 1, ivx) * sk(i, j, k - 1, 3)) * fact

        ! Idem for the gradient of v.

        dvdx = (w(i + 1, j, k, ivy) * si(i, j, k, 1) - w(i - 1, j, k, ivy) * si(i - 1, j, k, 1) &
              + w(i, j + 1, k, ivy) * sj(i, j, k, 1) - w(i, j - 1, k, ivy) * sj(i, j - 1, k, 1) &
              + w(i, j, k + 1, ivy) * sk(i, j, k, 1) - w(i, j, k - 1, ivy) * sk(i, j, k - 1, 1)) * fact
        dvdy = (w(i + 1, j, k, ivy) * si(i, j, k, 2) - w(i - 1, j, k, ivy) * si(i - 1, j, k, 2) &
              + w(i, j + 1, k, ivy) * sj(i, j, k, 2) - w(i, j - 1, k, ivy) * sj(i, j - 1, k, 2) &
              + w(i, j, k + 1, ivy) * sk(i, j, k, 2) - w(i, j, k - 1, ivy) * sk(i, j, k - 1, 2)) * fact
        dvdz = (w(i + 1, j, k, ivy) * si(i, j, k, 3) - w(i - 1, j, k, ivy) * si(i - 1, j, k, 3) &
              + w(i, j + 1, k, ivy) * sj(i, j, k, 3) - w(i, j - 1, k, ivy) * sj(i, j - 1, k, 3) &
              + w(i, j, k + 1, ivy) * sk(i, j, k, 3) - w(i, j, k - 1, ivy) * sk(i, j, k - 1, 3)) * fact

        ! And for the gradient of w.

        dwdx = (w(i + 1, j, k, ivz) * si(i, j, k, 1) - w(i - 1, j, k, ivz) * si(i - 1, j, k, 1) &
              + w(i, j + 1, k, ivz) * sj(i, j, k, 1) - w(i, j - 1, k, ivz) * sj(i, j - 1, k, 1) &
              + w(i, j, k + 1, ivz) * sk(i, j, k, 1) - w(i, j, k - 1, ivz) * sk(i, j, k - 1, 1)) * fact
        dwdy = (w(i + 1, j, k, ivz) * si(i, j, k, 2) - w(i - 1, j, k, ivz) * si(i - 1, j, k, 2) &
              + w(i, j + 1, k, ivz) * sj(i, j, k, 2) - w(i, j - 1, k, ivz) * sj(i, j - 1, k, 2) &
              + w(i, j, k + 1, ivz) * sk(i, j, k, 2) - w(i, j, k - 1, ivz) * sk(i, j, k - 1, 2)) * fact
        dwdz = (w(i + 1, j, k, ivz) * si(i, j, k, 3) - w(i - 1, j, k, ivz) * si(i - 1, j, k, 3) &
              + w(i, j + 1, k, ivz) * sj(i, j, k, 3) - w(i, j - 1, k, ivz) * sj(i, j - 1, k, 3) &
              + w(i, j, k + 1, ivz) * sk(i, j, k, 3) - w(i, j, k - 1, ivz) * sk(i, j, k - 1, 3)) * fact
        
        dU_dx = U_inv*(w(i, j, k, ivx)*dudx + w(i, j, k, ivy)*dudy + w(i, j, k, ivz)*dudz)
        dU_dy = U_inv*(w(i, j, k, ivx)*dvdx + w(i, j, k, ivy)*dvdy + w(i, j, k, ivz)*dvdz)
        dU_dz = U_inv*(w(i, j, k, ivx)*dwdx + w(i, j, k, ivy)*dwdy + w(i, j, k, ivz)*dwdz)

        dU_ds = w(i, j, k, ivx)/U*dU_dx + &
            w(i, j, k, ivy)/U*dU_dy + &
            w(i, j, k, ivz)/U*dU_dz

        Tu = 100.0 * sqrt(2*w(i, j, k, itu1)/3) / U
        Tu = max(Tu, 0.027) ! clip for numerical robustness


        ! Now we need to solve for theta through Newton's method. The number of iterations is hard-coded so tapenade is 
        ! able to differentiate it
        thetat = 0.01
        do n = 1, 10
            lambda = w(i, j, k, irho)*thetat**2 / rlv(i, j, k) * dU_ds
            lambda = max(min(lambda, 0.1), -0.1) ! clip for numerical robustness

            ! compute F function
            F1 = 1.0 + 0.275*(1.0 - exp(-35.0 * lambda))*exp(-(Tu/0.5))
            F2 = max(F1, 1.0)
            F3 = 1.0 - (-12.986*lambda - 123.66*lambda**2 - 405.689*lambda**3) * &
                exp(-(Tu/1.5)**1.5)
            F = min(F2, F3)

            if (Tu .gt. 1.3) then
                Re_thetat_eq_1 = 331.50*(Tu-0.5658)**(-0.671) * F
            else
                Re_thetat_eq_1 = (1173.51 - 589.428*Tu + 0.2196*Tu**(-2)) * F
            end if

            Re_thetat_eq_2 = w(i, j, k, irho) * U *thetat / rlv(i, j, k)

            ! residum which should go to 0
            residum = Re_thetat_eq_1 - Re_thetat_eq_2

            ! print *, 'iteration, thetat, residum, re_thetat_eq_1', n, thetat, residum, re_thetat_eq_1

            ! if we are in the first iteration, we need to kickstart the secant method first
            if (n .eq. 1) then
                residum_old = residum
                thetat_old = thetat
                thetat = 0.5*thetat
                cycle
            end if

            ! if the residum is basically 0, we cycle until we reach the end. (we cant drop out because tapenade would be unable to
            ! differentiate it)
            if ( abs(residum) .lt. 1e-9) then
                cycle
            end if

            ! compute next step (secant method)
            thetat_new = (thetat_old*residum - thetat*residum_old) / (residum - residum_old)

            ! save values for next iteration
            residum_old = residum
            thetat_old = thetat
            thetat = thetat_new

        end do

        ! save result in output-variable
        Re_thetat_eq = max(Re_thetat_eq_1, 20.0) ! clip for numerical robustness

        ! print *, 'thetat, lambda, Re_thetat_eq', thetat, lambda, Re_thetat_eq

    end subroutine solve_local_Re_thetat_eq


    subroutine GammaRethetaSource
        use blockPointers
        use constants
        use variableConstants
        use paramTurb
        implicit None

        integer(kind=intType) :: i, j, k
        real(kind=realType) :: Re_thetat_eq, U2, U, lambda_theta, delta, F_theta_t, T, R_t, Re_theta_c
        real(kind=realType) :: Re_S, F_length1, F_length, F_onset1, F_onset, F_turb, P_gamma, E_gamma, P_thetat
        real(kind=realType) :: Re_omega, F_wake

        real(kind=realType) :: rhoi, vort


#ifdef TAPENADE_REVERSE
        !$AD II-LOOP
        do ii = 0, nx * ny * nz - 1
            i = mod(ii, nx) + 2
            j = mod(ii / nx, ny) + 2
            k = ii / (nx * ny) + 2
#else
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il
#endif
                        rhoi = one / w(i, j, k, irho)

                        vort = max(sqrt(scratch(i, j, k, iVorticity)), eps)

                        ! compute Re_thetat_eq
                        call solve_local_Re_thetat_eq(Re_thetat_eq, i, j, k)
                        
                        U2 = w(i, j, k, ivx)**2 + w(i, j, k, ivy)**2 + w(i, j, k, ivz)**2
                        U = sqrt(U2)

                        Re_omega = (w(i, j, k, irho)*w(i, j, k, itu2)*d2wall(i, j, k)**2) / rlv(i, j, k)
                        F_wake = exp(-(Re_omega/1e5)**2)

                        delta = 375.0*vort*w(i, j, k, iTransition2)*d2wall(i, j, k) / &
                                (w(i, j, k, irho) * U)
                        F_theta_t = min(max(F_wake*exp(-(d2wall(i, j, k)/delta)**4), & ! todo: pull out of scratch
                                1.0 - ((rLMce2*w(i, j, k, iTransition1) - 1.0)/(rLMce2-1))**2), 1.0)


                        T = 500.0 * rlv(i, j, k) / (w(i, j, k, irho) * U2)

                        
                        R_t = w(i, j, k, irho) * w(i, j, k, itu1) / (rlv(i, j, k) * w(i, j, k, itu2)) ! todo: save this in scratch
                        Re_theta_c = 0.67*w(i, j, k, iTransition2) + 24.0*sin(w(i, j, k, iTransition2)/240.0 + 0.5) + 14.0 ! todo: save this in scratch
                        Re_S = w(i, j, k, irho) * sqrt(scratch(i, j, k, iStrain)) * d2wall(i, j, k)**2 / rev(i, j, k) ! todo: save this in scratch

                        F_length1 = exp(-0.03 * (w(i, j, k, iTransition2) - 460.0))
                        F_length = 44.0 - (44.0 - (0.5 - 30000.0 * (w(i, j, k, iTransition2) - 596.0))) / &
                            (1 + F_length1)**(1/6)

                        F_onset1 = sqrt(Re_S /(2.6 * Re_theta_c)) ! continue here
                        F_onset = (tanh(6.0*(F_onset1 - 1.35)) + 1.0) / 2.0
                        F_turb = (1 - F_onset) * exp(-R_t)


                        ! since we need to divide by rho, rho does not appear here anymore
                        P_gamma = F_length * rLMca1 * w(i, j, k, irho) * sqrt(scratch(i, j, k, iStrain)) * &
                            sqrt(w(i, j, k, iTransition1)* F_onset)*(1.0 - rLMce1*w(i, j, k, iTransition1))

                        E_gamma = rLMca2 * w(i, j, k, irho) * vort * w(i, j, k, iTransition1) * &
                            F_turb * (rLMce2 * w(i, j, k, iTransition1) - 1.0)
                            
                        P_thetat = rLMcthetat * w(i, j, k, irho) / T * (Re_thetat_eq - w(i, j, k, iTransition2)) * &
                            (1.0 - F_theta_t)

                        scratch(i, j, k, isTransition1) = (P_gamma - E_gamma) * rhoi
                        scratch(i, j, k, isTransition2) = P_thetat * rhoi

                        ! print *, 'source terms: gamma, thetat', scratch(i, j, k, isTransition1), scratch(i, j, k, isTransition2)

#ifdef TAPENADE_REVERSE
                    end do
#else
                end do
            end do
        end do
#endif
    end subroutine GammaRethetaSource

    subroutine GammaRethetaViscous

        use blockPointers
        use constants
        use variableConstants
        use paramTurb
        implicit none
        !
        !      Local variables.
        !
        integer(kind=intType) :: i, j, k, ii

        real(kind=realType) :: rhoi
        real(kind=realType) :: voli, volmi, volpi
        real(kind=realType) :: xm, ym, zm, xp, yp, zp, xa, ya, za
        real(kind=realType) :: mulm, mulp, muem, muep
        real(kind=realType) :: ttm, ttp
        real(kind=realType) :: c1m, c1p, c10, c2m, c2p, c20
        real(kind=realType) :: rblank

        !       Advection and unsteady terms.
        !
        !
        !       Viscous terms in k-direction.
        !
#ifdef TAPENADE_REVERSE
        !$AD II-LOOP
        do ii = 0, nx * ny * nz - 1
            i = mod(ii, nx) + 2
            j = mod(ii / nx, ny) + 2
            k = ii / (nx * ny) + 2
#else
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il
#endif

                        ! Compute the metrics in zeta-direction, i.e. along the
                        ! line k = constant.

                        voli = one / vol(i, j, k)
                        volmi = two / (vol(i, j, k) + vol(i, j, k - 1))
                        volpi = two / (vol(i, j, k) + vol(i, j, k + 1))

                        xm = sk(i, j, k - 1, 1) * volmi
                        ym = sk(i, j, k - 1, 2) * volmi
                        zm = sk(i, j, k - 1, 3) * volmi
                        xp = sk(i, j, k, 1) * volpi
                        yp = sk(i, j, k, 2) * volpi
                        zp = sk(i, j, k, 3) * volpi

                        xa = half * (sk(i, j, k, 1) + sk(i, j, k - 1, 1)) * voli
                        ya = half * (sk(i, j, k, 2) + sk(i, j, k - 1, 2)) * voli
                        za = half * (sk(i, j, k, 3) + sk(i, j, k - 1, 3)) * voli
                        ttm = xm * xa + ym * ya + zm * za
                        ttp = xp * xa + yp * ya + zp * za

                        ! Computation of the viscous terms in zeta-direction; note
                        ! that cross-derivatives are neglected, i.e. the mesh is
                        ! assumed to be orthogonal.
                        ! The second derivative in zeta-direction is constructed as
                        ! the central difference of the first order derivatives, i.e.
                        ! d^2/dzeta^2 = d/dzeta (d/dzeta k+1/2 - d/dzeta k-1/2).
                        ! In this way the metric as well as the varying viscosity
                        ! can be taken into account; the latter appears inside the
                        ! d/dzeta derivative. The whole term is divided by rho to
                        ! obtain the diffusion term for k and omega.

                        ! First the gamma-term.

                        rhoi = one / w(i, j, k, irho)
                        mulm = half * (rlv(i, j, k - 1) + rlv(i, j, k))
                        mulp = half * (rlv(i, j, k + 1) + rlv(i, j, k))
                        muem = half * (rev(i, j, k - 1) + rev(i, j, k))
                        muep = half * (rev(i, j, k + 1) + rev(i, j, k))

                        c1m = ttm * (mulm + muem / rLMsigmaf) * rhoi
                        c1p = ttp * (mulp + muep / rLMsigmaf) * rhoi
                        c10 = c1m + c1p

                        ! And the re_theta_t term.

                        muem = half * (rev(i, j, k - 1) + rev(i, j, k))
                        muep = half * (rev(i, j, k + 1) + rev(i, j, k))

                        c2m = ttm * rLMsigmathetat * (mulm + muem) * rhoi
                        c2p = ttp * rLMsigmathetat * (mulp + muep) * rhoi
                        c20 = c2m + c2p

                        ! Update the residual for this cell and store the possible
                        ! coefficients for the matrix in b1, b2, c1, c2, d1 and d2.

                        scratch(i, j, k, isTransition1) = scratch(i, j, k, isTransition1) + &
                            c1m * w(i, j, k - 1, iTransition1) - &
                            c10 * w(i, j, k, iTransition1) + &
                            c1p * w(i, j, k + 1, iTransition1)

                        scratch(i, j, k, isTransition2) = scratch(i, j, k, isTransition2) + &
                            c2m * w(i, j, k - 1, iTransition2) - &
                            c20 * w(i, j, k, iTransition2) + &
                            c2p * w(i, j, k + 1, iTransition2)
! #ifndef USE_TAPENADE
!                         b1 = -c1m
!                         c1 = c10
!                         d1 = -c1p

!                         b2 = -c2m
!                         c2 = c20
!                         d2 = -c2p

!                         ! Update the central jacobian. For nonboundary cells this
!                         ! is simply c1 and c2. For boundary cells this is slightly
!                         ! more complicated, because the boundary conditions are
!                         ! treated implicitly and the off-diagonal terms b1, b2 and
!                         ! d1, d2 must be taken into account.
!                         ! The boundary conditions are only treated implicitly if
!                         ! the diagonal dominance of the matrix is increased.

!                         if (k == 2) then
!                             qq(i, j, k, 1, 1) = qq(i, j, k, 1, 1) + c1 &
!                                                 - b1 * max(bmtk1(i, j, itu1, itu1), zero)
!                             qq(i, j, k, 1, 2) = qq(i, j, k, 1, 2) - b1 * bmtk1(i, j, itu1, itu2)
!                             qq(i, j, k, 2, 1) = qq(i, j, k, 2, 1) - b2 * bmtk1(i, j, itu2, itu1)
!                             qq(i, j, k, 2, 2) = qq(i, j, k, 2, 2) + c2 &
!                                                 - b2 * max(bmtk1(i, j, itu2, itu2), zero)
!                         else if (k == kl) then
!                             qq(i, j, k, 1, 1) = qq(i, j, k, 1, 1) + c1 &
!                                                 - d1 * max(bmtk2(i, j, itu1, itu1), zero)
!                             qq(i, j, k, 1, 2) = qq(i, j, k, 1, 2) - d1 * bmtk2(i, j, itu1, itu2)
!                             qq(i, j, k, 2, 1) = qq(i, j, k, 2, 1) - d2 * bmtk2(i, j, itu2, itu1)
!                             qq(i, j, k, 2, 2) = qq(i, j, k, 2, 2) + c2 &
!                                                 - d2 * max(bmtk2(i, j, itu2, itu2), zero)
!                         else
!                             qq(i, j, k, 1, 1) = qq(i, j, k, 1, 1) + c1
!                             qq(i, j, k, 2, 2) = qq(i, j, k, 2, 2) + c2
!                         end if
! #endif
#ifdef TAPENADE_REVERSE
                    end do
#else
                end do
            end do
        end do
#endif
        !
        !       Viscous terms in j-direction.
        !
#ifdef TAPENADE_REVERSE
        !$AD II-LOOP
        do ii = 0, nx * ny * nz - 1
            i = mod(ii, nx) + 2
            j = mod(ii / nx, ny) + 2
            k = ii / (nx * ny) + 2
#else
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il
#endif

                        ! Compute the metrics in eta-direction, i.e. along the
                        ! line j = constant.

                        voli = one / vol(i, j, k)
                        volmi = two / (vol(i, j, k) + vol(i, j - 1, k))
                        volpi = two / (vol(i, j, k) + vol(i, j + 1, k))

                        xm = sj(i, j - 1, k, 1) * volmi
                        ym = sj(i, j - 1, k, 2) * volmi
                        zm = sj(i, j - 1, k, 3) * volmi
                        xp = sj(i, j, k, 1) * volpi
                        yp = sj(i, j, k, 2) * volpi
                        zp = sj(i, j, k, 3) * volpi

                        xa = half * (sj(i, j, k, 1) + sj(i, j - 1, k, 1)) * voli
                        ya = half * (sj(i, j, k, 2) + sj(i, j - 1, k, 2)) * voli
                        za = half * (sj(i, j, k, 3) + sj(i, j - 1, k, 3)) * voli
                        ttm = xm * xa + ym * ya + zm * za
                        ttp = xp * xa + yp * ya + zp * za

                        ! Computation of the viscous terms in eta-direction; note
                        ! that cross-derivatives are neglected, i.e. the mesh is
                        ! assumed to be orthogonal.
                        ! The second derivative in eta-direction is constructed as
                        ! the central difference of the first order derivatives, i.e.
                        ! d^2/deta^2 = d/deta (d/deta j+1/2 - d/deta j-1/2).
                        ! In this way the metric as well as the varying viscosity
                        ! can be taken into account; the latter appears inside the
                        ! d/deta derivative. The whole term is divided by rho to
                        ! obtain the diffusion term for k and omega.

                        ! First the gamma-term.

                        rhoi = one / w(i, j, k, irho)
                        mulm = half * (rlv(i, j - 1, k) + rlv(i, j, k))
                        mulp = half * (rlv(i, j + 1, k) + rlv(i, j, k))
                        muem = half * (rev(i, j - 1, k) + rev(i, j, k))
                        muep = half * (rev(i, j + 1, k) + rev(i, j, k))

                        c1m = ttm * (mulm + muem / rLMsigmaf) * rhoi
                        c1p = ttp * (mulp + muep / rLMsigmaf) * rhoi
                        c10 = c1m + c1p

                        ! And the re_theta_t term.

                        muem = half * (rev(i, j - 1, k) + rev(i, j, k))
                        muep = half * (rev(i, j + 1, k) + rev(i, j, k))

                        c2m = ttm * rLMsigmathetat * (mulm + muem) * rhoi
                        c2p = ttp * rLMsigmathetat * (mulp + muep) * rhoi
                        c20 = c2m + c2p

                        ! Update the residual for this cell and store the possible
                        ! coefficients for the matrix in b1, b2, c1, c2, d1 and d2.

                        scratch(i, j, k, isTransition1) = scratch(i, j, k, isTransition1) + &
                            c1m * w(i, j - 1, k, iTransition1) - &
                            c10 * w(i, j, k, iTransition1) + &
                            c1p * w(i, j + 1, k, iTransition1)

                        scratch(i, j, k, isTransition2) = scratch(i, j, k, isTransition2) + &
                            c2m * w(i, j - 1, k, iTransition2) - &
                            c20 * w(i, j, k, iTransition2) + &
                            c2p * w(i, j + 1, k, iTransition2)

! #ifndef USE_TAPENADE
!                         b1 = -c1m
!                         c1 = c10
!                         d1 = -c1p

!                         b2 = -c2m
!                         c2 = c20
!                         d2 = -c2p

!                         ! Update the central jacobian. For nonboundary cells this
!                         ! is simply c1 and c2. For boundary cells this is slightly
!                         ! more complicated, because the boundary conditions are
!                         ! treated implicitly and the off-diagonal terms b1, b2 and
!                         ! d1, d2 must be taken into account.
!                         ! The boundary conditions are only treated implicitly if
!                         ! the diagonal dominance of the matrix is increased.

!                         if (j == 2) then
!                             qq(i, j, k, 1, 1) = qq(i, j, k, 1, 1) + c1 &
!                                                 - b1 * max(bmtj1(i, k, itu1, itu1), zero)
!                             qq(i, j, k, 1, 2) = qq(i, j, k, 1, 2) - b1 * bmtj1(i, k, itu1, itu2)
!                             qq(i, j, k, 2, 1) = qq(i, j, k, 2, 1) - b2 * bmtj1(i, k, itu2, itu1)
!                             qq(i, j, k, 2, 2) = qq(i, j, k, 2, 2) + c2 &
!                                                 - b2 * max(bmtj1(i, k, itu2, itu2), zero)
!                         else if (j == jl) then
!                             qq(i, j, k, 1, 1) = qq(i, j, k, 1, 1) + c1 &
!                                                 - d1 * max(bmtj2(i, k, itu1, itu1), zero)
!                             qq(i, j, k, 1, 2) = qq(i, j, k, 1, 2) - d1 * bmtj2(i, k, itu1, itu2)
!                             qq(i, j, k, 2, 1) = qq(i, j, k, 2, 1) - d2 * bmtj2(i, k, itu2, itu1)
!                             qq(i, j, k, 2, 2) = qq(i, j, k, 2, 2) + c2 &
!                                                 - d2 * max(bmtj2(i, k, itu2, itu2), zero)
!                         else
!                             qq(i, j, k, 1, 1) = qq(i, j, k, 1, 1) + c1
!                             qq(i, j, k, 2, 2) = qq(i, j, k, 2, 2) + c2
!                         end if

! #endif
#ifdef TAPENADE_REVERSE
                    end do
#else
                end do
            end do
        end do
#endif
        !
        !       Viscous terms in i-direction.
        !
#ifdef TAPENADE_REVERSE
        !$AD II-LOOP
        do ii = 0, nx * ny * nz - 1
            i = mod(ii, nx) + 2
            j = mod(ii / nx, ny) + 2
            k = ii / (nx * ny) + 2
#else
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il
#endif
                        ! Compute the metrics in xi-direction, i.e. along the
                        ! line i = constant.

                        voli = one / vol(i, j, k)
                        volmi = two / (vol(i, j, k) + vol(i - 1, j, k))
                        volpi = two / (vol(i, j, k) + vol(i + 1, j, k))

                        xm = si(i - 1, j, k, 1) * volmi
                        ym = si(i - 1, j, k, 2) * volmi
                        zm = si(i - 1, j, k, 3) * volmi
                        xp = si(i, j, k, 1) * volpi
                        yp = si(i, j, k, 2) * volpi
                        zp = si(i, j, k, 3) * volpi

                        xa = half * (si(i, j, k, 1) + si(i - 1, j, k, 1)) * voli
                        ya = half * (si(i, j, k, 2) + si(i - 1, j, k, 2)) * voli
                        za = half * (si(i, j, k, 3) + si(i - 1, j, k, 3)) * voli
                        ttm = xm * xa + ym * ya + zm * za
                        ttp = xp * xa + yp * ya + zp * za

                        ! Computation of the viscous terms in xi-direction; note
                        ! that cross-derivatives are neglected, i.e. the mesh is
                        ! assumed to be orthogonal.
                        ! The second derivative in xi-direction is constructed as
                        ! the central difference of the first order derivatives, i.e.
                        ! d^2/dxi^2 = d/dxi (d/dxi i+1/2 - d/dxi i-1/2).
                        ! In this way the metric as well as the varying viscosity
                        ! can be taken into account; the latter appears inside the
                        ! d/dxi derivative. The whole term is divided by rho to
                        ! obtain the diffusion term for k and omega.

                        ! First the gamma-term.

                        rhoi = one / w(i, j, k, irho)
                        mulm = half * (rlv(i - 1, j, k) + rlv(i, j, k))
                        mulp = half * (rlv(i + 1, j, k) + rlv(i, j, k))
                        muem = half * (rev(i - 1, j, k) + rev(i, j, k))
                        muep = half * (rev(i + 1, j, k) + rev(i, j, k))

                        c1m = ttm * (mulm + muem / rLMsigmaf) * rhoi
                        c1p = ttp * (mulp + muep / rLMsigmaf) * rhoi
                        c10 = c1m + c1p

                        ! And the re_theta_t term.

                        muem = half * (rev(i - 1, j, k) + rev(i, j, k))
                        muep = half * (rev(i + 1, j, k) + rev(i, j, k))

                        c2m = ttm * rLMsigmathetat * (mulm + muem) * rhoi
                        c2p = ttp * rLMsigmathetat * (mulp + muep) * rhoi
                        c20 = c2m + c2p

                        ! Update the residual for this cell and store the possible
                        ! coefficients for the matrix in b1, b2, c1, c2, d1 and d2.

                        scratch(i, j, k, isTransition1) = scratch(i, j, k, isTransition1) + &
                            c1m * w(i - 1, j, k, iTransition1) - &
                            c10 * w(i, j, k, iTransition1) + &
                            c1p * w(i + 1, j, k, iTransition1)

                        scratch(i, j, k, isTransition2) = scratch(i, j, k, isTransition2) + &
                            c2m * w(i - 1, j, k, iTransition2) - &
                            c20 * w(i, j, k, iTransition2) + &
                            c2p * w(i + 1, j, k, iTransition2)

! #ifndef USE_TAPENADE
!                         b1 = -c1m
!                         c1 = c10
!                         d1 = -c1p

!                         b2 = -c2m
!                         c2 = c20
!                         d2 = -c2p

!                         ! Update the central jacobian. For nonboundary cells this
!                         ! is simply c1 and c2. For boundary cells this is slightly
!                         ! more complicated, because the boundary conditions are
!                         ! treated implicitly and the off-diagonal terms b1, b2 and
!                         ! d1, d2 must be taken into account.
!                         ! The boundary conditions are only treated implicitly if
!                         ! the diagonal dominance of the matrix is increased.

!                         if (i == 2) then
!                             qq(i, j, k, 1, 1) = qq(i, j, k, 1, 1) + c1 &
!                                                 - b1 * max(bmti1(j, k, itu1, itu1), zero)
!                             qq(i, j, k, 1, 2) = qq(i, j, k, 1, 2) - b1 * bmti1(j, k, itu1, itu2)
!                             qq(i, j, k, 2, 1) = qq(i, j, k, 2, 1) - b2 * bmti1(j, k, itu2, itu1)
!                             qq(i, j, k, 2, 2) = qq(i, j, k, 2, 2) + c2 &
!                                                 - b2 * max(bmti1(j, k, itu2, itu2), zero)
!                         else if (i == il) then
!                             qq(i, j, k, 1, 1) = qq(i, j, k, 1, 1) + c1 &
!                                                 - d1 * max(bmti2(j, k, itu1, itu1), zero)
!                             qq(i, j, k, 1, 2) = qq(i, j, k, 1, 2) - d1 * bmti2(j, k, itu1, itu2)
!                             qq(i, j, k, 2, 1) = qq(i, j, k, 2, 1) - d2 * bmti2(j, k, itu2, itu1)
!                             qq(i, j, k, 2, 2) = qq(i, j, k, 2, 2) + c2 &
!                                                 - d2 * max(bmti2(j, k, itu2, itu2), zero)
!                         else
!                             qq(i, j, k, 1, 1) = qq(i, j, k, 1, 1) + c1
!                             qq(i, j, k, 2, 2) = qq(i, j, k, 2, 2) + c2
!                         end if
! #endif
#ifdef TAPENADE_REVERSE
                    end do
#else
                end do
            end do
        end do
#endif


    end subroutine GammaRethetaViscous

    subroutine GammaRethetaResScale
        use blockPointers
        use constants
        use variableConstants

        implicit none
        !
        !      Local variables.
        !
        integer(kind=intType) :: i, j, k, ii

        real(kind=realType) :: rblank

        ! Multiply the residual by the volume and store this in dw; this
        ! is done for monitoring reasons only. The multiplication with the
        ! volume is present to be consistent with the flow residuals; also
        ! the negative value is taken, again to be consistent with the
        ! flow equations. Also multiply by iblank so that no updates occur
        ! in holes or the overset boundary.

#ifdef TAPENADE_REVERSE
        !$AD II-LOOP
        do ii = 0, nx * ny * nz - 1
            i = mod(ii, nx) + 2
            j = mod(ii / nx, ny) + 2
            k = ii / (nx * ny) + 2
#else
            do k = 2, kl
                do j = 2, jl
                    do i = 2, il
#endif
                        rblank = real(iblank(i, j, k), realType)
                        dw(i, j, k, iTransition1) = -volRef(i, j, k) * scratch(i, j, k, isTransition1) * rblank
                        dw(i, j, k, iTransition2) = -volRef(i, j, k) * scratch(i, j, k, isTransition2) * rblank

                        ! print *, 'dw Transition 1 & 2', dw(i, j, k, iTransition1), dw(i, j, k, iTransition2)
#ifdef TAPENADE_REVERSE
                    end do
#else
                end do
            end do
        end do
#endif

    end subroutine GammaRethetaResScale


end module GammaRethetaModel
