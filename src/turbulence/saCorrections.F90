! This module contains the source terms for two variants of the
! Spalart-Allmaras turbulence model that are not part of the
! "standard" SA implementation in sa.F90:
!
!   1. The Edwards-Chandra modification (SA-noft2-Edwards), selected with
!      turbulenceModel = "SA-Edwards":
!
!        Edwards, J. R. and Chandra, S., "Comparison of Eddy
!        Viscosity-Transport Turbulence Models for Three-Dimensional,
!        Shock-Separated Flowfields," AIAA Paper 94-2275, 1994,
!        https://doi.org/10.2514/6.1994-2275; also AIAA Journal,
!        Vol. 34, No. 4, 1996, pp. 756-763, https://doi.org/10.2514/3.13137.
!
!      The model is identical to standard SA except that ft2 is ignored
!      and
!
!        Stilde = sqrt(S) * (1/chi + fv1)
!        r      = tanh(nuTilde / (Stilde kappa^2 d^2)) / tanh(1)
!        S      = (du_i/dx_j + du_j/dx_i) du_i/dx_j - 2/3 (du_k/dx_k)^2
!
!      i.e. the strain rate, not the vorticity, drives production. Note
!      that Stilde * nuTilde = sqrt(S) * (nu + nuTilde * fv1), which is how
!      the terms are evaluated here to avoid the 1/chi singularity at walls.
!
!   2. The mixing-layer compressibility correction (SA-comp), selected with
!      useCompressibilitySA = True, which can be combined with either
!      "SA" or "SA-Edwards":
!
!        Spalart, P. R., "Trends in Turbulence Treatments," AIAA Paper
!        2000-2306, 2000, https://doi.org/10.2514/6.2000-2306.
!
!      The following term is added to the right hand side of the
!      (non-conservative) nuTilde equation,
!
!        -C5 * nuTilde^2 / a^2 * (du_i/dx_j) (du_i/dx_j),   C5 = 3.5,
!
!      where a is the local speed of sound.
!
! The routines follow the conventions of saSource in sa.F90: the source
! is stored in scratch(:,:,:,idvt) and -dSource/dnuTilde is stored in qq
! for the DDADI solver. They are differentiated with Tapenade (see
! adjoint/Makefile_tapenade).

module saCorrections

    use constants
    implicit none

contains
#ifndef USE_TAPENADE
    subroutine saCorr_block(resOnly)
        !
        !       saCorr_block is the counterpart of sa_block for the SA
        !       variants in this module. It uses the Edwards source term
        !       for turbModel == spalartAllmarasEdwards and the standard
        !       source term otherwise, and adds the compressibility
        !       correction if useCompressibilitySA is set.
        !
        use constants
        use blockPointers, only: il, jl, kl
        use inputPhysics, only: turbModel, useCompressibilitySA
        use sa, only: qq, saSource, saViscous, saResScale, saSolve
        use turbutils, only: turbAdvection, unsteadyTurbTerm, saEddyViscosity
        use turbBCRoutines, only: bcTurbTreatment, applyAllTurbBCThisBlock
        implicit none
        !
        !      Subroutine argument.
        !
        logical, intent(in) :: resOnly
        !
        !      Local variables.
        !
        integer(kind=intType) :: nn

        ! Set the arrays for the boundary condition treatment.
        call bcTurbTreatment

        ! Alloc central jacobian memory
        allocate (qq(2:il, 2:jl, 2:kl))

        ! Source Terms
        if (turbModel == spalartAllmarasEdwards) then
            call saEdwardsSource
        else
            call saSource
        end if

        if (useCompressibilitySA) then
            call saCompressibilitySource
        end if

        ! Advection Term
        nn = itu1 - 1
        call turbAdvection(1_intType, 1_intType, nn, qq)

        ! Unsteady Term
        call unsteadyTurbTerm(1_intType, 1_intType, nn, qq)

        ! Viscous Terms
        call saViscous

        ! Perform the residual scaling
        call saResScale

        if (.not. resOnly) then

            ! Do solve
            call saSolve

            ! Compute the corresponding eddy viscosity.
            call saEddyViscosity(2, il, 2, jl, 2, kl)

            ! Set the halo values for the turbulent variables.
            call applyAllTurbBCThisBlock(.true.)
        end if

        deallocate (qq)
    end subroutine saCorr_block
#endif

    subroutine saEdwardsSource
        !
        !  Source terms of the Edwards-Chandra SA model.
        !  Determine the source term and its derivative w.r.t. nuTilde
        !  for all internal cells of the block.
        !  Remember that the SA field variable nuTilde = w(i,j,k,itu1)

        use blockPointers
        use constants
        use paramTurb
        use inputDiscretization, only: approxSA
#ifndef USE_TAPENADE
        use sa, only: qq
#endif
        implicit none

        ! Local parameters
        real(kind=realType), parameter :: f23 = two * third
        real(kind=realType), parameter :: xminn = 1.e-10_realType

        ! Local variables.
        integer(kind=intType) :: i, j, k, ii
        real(kind=realType) :: cv13, kar2Inv, cw36, tanhOneInv
        real(kind=realType) :: fv1, nuFv, ss, nu, dist2Inv, chi, chi2, chi3
        real(kind=realType) :: xx, tanhX, rr, gg, gg6, termFw, fwSa
        real(kind=realType) :: prod, dest
        real(kind=realType) :: dfv1, dnuFv, dxx, drr, dgg, dfw
        real(kind=realType) :: uux, uuy, uuz, vvx, vvy, vvz, wwx, wwy, wwz
        real(kind=realType) :: div2, fact, sxx, syy, szz, sxy, sxz, syz
        real(kind=realType) :: strainMag2, strainProd

        ! Set model constants
        cv13 = rsaCv1**3
        kar2Inv = one / (rsaK**2)
        cw36 = rsaCw3**6
        tanhOneInv = one / tanh(one)

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
                        ! Compute the velocity gradients in the cell center,
                        ! scaled by 2*vol. See saSource in sa.F90.

                        uux = w(i + 1, j, k, ivx) * si(i, j, k, 1) - w(i - 1, j, k, ivx) * si(i - 1, j, k, 1) &
                              + w(i, j + 1, k, ivx) * sj(i, j, k, 1) - w(i, j - 1, k, ivx) * sj(i, j - 1, k, 1) &
                              + w(i, j, k + 1, ivx) * sk(i, j, k, 1) - w(i, j, k - 1, ivx) * sk(i, j, k - 1, 1)
                        uuy = w(i + 1, j, k, ivx) * si(i, j, k, 2) - w(i - 1, j, k, ivx) * si(i - 1, j, k, 2) &
                              + w(i, j + 1, k, ivx) * sj(i, j, k, 2) - w(i, j - 1, k, ivx) * sj(i, j - 1, k, 2) &
                              + w(i, j, k + 1, ivx) * sk(i, j, k, 2) - w(i, j, k - 1, ivx) * sk(i, j, k - 1, 2)
                        uuz = w(i + 1, j, k, ivx) * si(i, j, k, 3) - w(i - 1, j, k, ivx) * si(i - 1, j, k, 3) &
                              + w(i, j + 1, k, ivx) * sj(i, j, k, 3) - w(i, j - 1, k, ivx) * sj(i, j - 1, k, 3) &
                              + w(i, j, k + 1, ivx) * sk(i, j, k, 3) - w(i, j, k - 1, ivx) * sk(i, j, k - 1, 3)

                        vvx = w(i + 1, j, k, ivy) * si(i, j, k, 1) - w(i - 1, j, k, ivy) * si(i - 1, j, k, 1) &
                              + w(i, j + 1, k, ivy) * sj(i, j, k, 1) - w(i, j - 1, k, ivy) * sj(i, j - 1, k, 1) &
                              + w(i, j, k + 1, ivy) * sk(i, j, k, 1) - w(i, j, k - 1, ivy) * sk(i, j, k - 1, 1)
                        vvy = w(i + 1, j, k, ivy) * si(i, j, k, 2) - w(i - 1, j, k, ivy) * si(i - 1, j, k, 2) &
                              + w(i, j + 1, k, ivy) * sj(i, j, k, 2) - w(i, j - 1, k, ivy) * sj(i, j - 1, k, 2) &
                              + w(i, j, k + 1, ivy) * sk(i, j, k, 2) - w(i, j, k - 1, ivy) * sk(i, j, k - 1, 2)
                        vvz = w(i + 1, j, k, ivy) * si(i, j, k, 3) - w(i - 1, j, k, ivy) * si(i - 1, j, k, 3) &
                              + w(i, j + 1, k, ivy) * sj(i, j, k, 3) - w(i, j - 1, k, ivy) * sj(i, j - 1, k, 3) &
                              + w(i, j, k + 1, ivy) * sk(i, j, k, 3) - w(i, j, k - 1, ivy) * sk(i, j, k - 1, 3)

                        wwx = w(i + 1, j, k, ivz) * si(i, j, k, 1) - w(i - 1, j, k, ivz) * si(i - 1, j, k, 1) &
                              + w(i, j + 1, k, ivz) * sj(i, j, k, 1) - w(i, j - 1, k, ivz) * sj(i, j - 1, k, 1) &
                              + w(i, j, k + 1, ivz) * sk(i, j, k, 1) - w(i, j, k - 1, ivz) * sk(i, j, k - 1, 1)
                        wwy = w(i + 1, j, k, ivz) * si(i, j, k, 2) - w(i - 1, j, k, ivz) * si(i - 1, j, k, 2) &
                              + w(i, j + 1, k, ivz) * sj(i, j, k, 2) - w(i, j - 1, k, ivz) * sj(i, j - 1, k, 2) &
                              + w(i, j, k + 1, ivz) * sk(i, j, k, 2) - w(i, j, k - 1, ivz) * sk(i, j, k - 1, 2)
                        wwz = w(i + 1, j, k, ivz) * si(i, j, k, 3) - w(i - 1, j, k, ivz) * si(i - 1, j, k, 3) &
                              + w(i, j + 1, k, ivz) * sj(i, j, k, 3) - w(i, j - 1, k, ivz) * sj(i, j - 1, k, 3) &
                              + w(i, j, k + 1, ivz) * sk(i, j, k, 3) - w(i, j, k - 1, ivz) * sk(i, j, k - 1, 3)

                        ! Strain rate tensor; the factor 1/(4*vol) accounts for
                        ! the 2*vol scaling of the gradients.

                        fact = fourth / vol(i, j, k)

                        sxx = two * fact * uux
                        syy = two * fact * vvy
                        szz = two * fact * wwz

                        sxy = fact * (uuy + vvx)
                        sxz = fact * (uuz + wwx)
                        syz = fact * (vvz + wwy)

                        div2 = f23 * (sxx + syy + szz)**2

                        strainMag2 = two * (sxy**2 + sxz**2 + syz**2) &
                                     + sxx**2 + syy**2 + szz**2

                        ! S = 2 s_ij s_ij - 2/3 (div u)^2. Edwards uses sqrt(S).
                        ! Limit it from below (like sst in sa.F90), which also
                        ! keeps the derivative of the square root finite.

                        strainProd = two * strainMag2 - div2
                        ss = sqrt(max(strainProd, xminn**2))

                        ! Laminar kinematic viscosity, inverse wall distance
                        ! squared, chi and fv1.

                        nu = rlv(i, j, k) / w(i, j, k, irho)
                        dist2Inv = one / (d2Wall(i, j, k)**2)
                        chi = w(i, j, k, itu1) / nu
                        chi2 = chi * chi
                        chi3 = chi * chi2
                        fv1 = chi3 / (chi3 + cv13)

                        ! nuFv = nu + nuTilde*fv1, such that
                        ! Stilde * nuTilde = ss * nuFv.

                        nuFv = nu + w(i, j, k, itu1) * fv1

                        ! nuFv >= nu for nuTilde >= 0. Guard against the zero
                        ! crossing that can occur for (unphysical) large
                        ! negative nuTilde.

                        nuFv = max(nuFv, xminn * nu)

                        ! r = tanh(nuTilde / (Stilde kappa^2 d^2)) / tanh(1)
                        !   = tanh(nuTilde^2 / (ss nuFv kappa^2 d^2)) / tanh(1).
                        ! No cut-off is needed since r <= 1/tanh(1).

                        xx = w(i, j, k, itu1)**2 * kar2Inv * dist2Inv / (ss * nuFv)
                        tanhX = tanh(xx)
                        rr = tanhX * tanhOneInv
                        gg = rr + rsaCw2 * (rr**6 - rr)
                        gg6 = gg**6
                        termFw = ((one + cw36) / (gg6 + cw36))**sixth
                        fwSa = gg * termFw

                        ! Production cb1*Stilde*nuTilde and destruction.

                        if (approxSA) then
                            prod = zero
                        else
                            prod = rsaCb1 * ss * nuFv
                        end if
                        dest = rsaCw1 * fwSa * dist2Inv * w(i, j, k, itu1)**2

                        scratch(i, j, k, idvt) = prod - dest

#ifndef USE_TAPENADE
                        ! Jacobian of the source term w.r.t. nuTilde. As in
                        ! saSource, the production is treated explicitly to
                        ! preserve diagonal dominance; only the destruction
                        ! is linearized. Note that -dsource/dnu is stored.

                        dfv1 = three * chi2 * cv13 / ((chi3 + cv13)**2)
                        dnuFv = fv1 + chi * dfv1

                        dxx = kar2Inv * dist2Inv / (ss * nuFv) &
                              * (two * w(i, j, k, itu1) - w(i, j, k, itu1)**2 * dnuFv / nuFv)
                        drr = (one - tanhX**2) * tanhOneInv * dxx
                        dgg = (one - rsaCw2 + six * rsaCw2 * (rr**5)) * drr
                        dfw = (cw36 / (gg6 + cw36)) * termFw * dgg

                        qq(i, j, k) = rsaCw1 * dist2Inv * w(i, j, k, itu1) &
                                      * (two * fwSa + w(i, j, k, itu1) * dfw)

                        qq(i, j, k) = max(qq(i, j, k), zero)
#endif
#ifdef TAPENADE_REVERSE
                    end do
#else
                end do
            end do
        end do
#endif
    end subroutine saEdwardsSource

    subroutine saCompressibilitySource
        !
        !  Compressibility correction of Spalart (AIAA 2000-2306).
        !  The term -C5 nuTilde^2/a^2 (du_i/dx_j)(du_i/dx_j) is added to
        !  the source term stored in scratch(:,:,:,idvt) by saSource or
        !  saEdwardsSource, and its (positive) contribution to
        !  -dsource/dnuTilde is added to qq.

        use blockPointers
        use constants
        use inputPhysics, only: SAc5
#ifndef USE_TAPENADE
        use sa, only: qq
#endif
        implicit none

        ! Local variables.
        integer(kind=intType) :: i, j, k, ii
        real(kind=realType) :: uux, uuy, uuz, vvx, vvy, vvz, wwx, wwy, wwz
        real(kind=realType) :: fact, gradU2, a2Inv, compCoef

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
                        ! Velocity gradients in the cell center, scaled by
                        ! 2*vol. See saSource in sa.F90.

                        uux = w(i + 1, j, k, ivx) * si(i, j, k, 1) - w(i - 1, j, k, ivx) * si(i - 1, j, k, 1) &
                              + w(i, j + 1, k, ivx) * sj(i, j, k, 1) - w(i, j - 1, k, ivx) * sj(i, j - 1, k, 1) &
                              + w(i, j, k + 1, ivx) * sk(i, j, k, 1) - w(i, j, k - 1, ivx) * sk(i, j, k - 1, 1)
                        uuy = w(i + 1, j, k, ivx) * si(i, j, k, 2) - w(i - 1, j, k, ivx) * si(i - 1, j, k, 2) &
                              + w(i, j + 1, k, ivx) * sj(i, j, k, 2) - w(i, j - 1, k, ivx) * sj(i, j - 1, k, 2) &
                              + w(i, j, k + 1, ivx) * sk(i, j, k, 2) - w(i, j, k - 1, ivx) * sk(i, j, k - 1, 2)
                        uuz = w(i + 1, j, k, ivx) * si(i, j, k, 3) - w(i - 1, j, k, ivx) * si(i - 1, j, k, 3) &
                              + w(i, j + 1, k, ivx) * sj(i, j, k, 3) - w(i, j - 1, k, ivx) * sj(i, j - 1, k, 3) &
                              + w(i, j, k + 1, ivx) * sk(i, j, k, 3) - w(i, j, k - 1, ivx) * sk(i, j, k - 1, 3)

                        vvx = w(i + 1, j, k, ivy) * si(i, j, k, 1) - w(i - 1, j, k, ivy) * si(i - 1, j, k, 1) &
                              + w(i, j + 1, k, ivy) * sj(i, j, k, 1) - w(i, j - 1, k, ivy) * sj(i, j - 1, k, 1) &
                              + w(i, j, k + 1, ivy) * sk(i, j, k, 1) - w(i, j, k - 1, ivy) * sk(i, j, k - 1, 1)
                        vvy = w(i + 1, j, k, ivy) * si(i, j, k, 2) - w(i - 1, j, k, ivy) * si(i - 1, j, k, 2) &
                              + w(i, j + 1, k, ivy) * sj(i, j, k, 2) - w(i, j - 1, k, ivy) * sj(i, j - 1, k, 2) &
                              + w(i, j, k + 1, ivy) * sk(i, j, k, 2) - w(i, j, k - 1, ivy) * sk(i, j, k - 1, 2)
                        vvz = w(i + 1, j, k, ivy) * si(i, j, k, 3) - w(i - 1, j, k, ivy) * si(i - 1, j, k, 3) &
                              + w(i, j + 1, k, ivy) * sj(i, j, k, 3) - w(i, j - 1, k, ivy) * sj(i, j - 1, k, 3) &
                              + w(i, j, k + 1, ivy) * sk(i, j, k, 3) - w(i, j, k - 1, ivy) * sk(i, j, k - 1, 3)

                        wwx = w(i + 1, j, k, ivz) * si(i, j, k, 1) - w(i - 1, j, k, ivz) * si(i - 1, j, k, 1) &
                              + w(i, j + 1, k, ivz) * sj(i, j, k, 1) - w(i, j - 1, k, ivz) * sj(i, j - 1, k, 1) &
                              + w(i, j, k + 1, ivz) * sk(i, j, k, 1) - w(i, j, k - 1, ivz) * sk(i, j, k - 1, 1)
                        wwy = w(i + 1, j, k, ivz) * si(i, j, k, 2) - w(i - 1, j, k, ivz) * si(i - 1, j, k, 2) &
                              + w(i, j + 1, k, ivz) * sj(i, j, k, 2) - w(i, j - 1, k, ivz) * sj(i, j - 1, k, 2) &
                              + w(i, j, k + 1, ivz) * sk(i, j, k, 2) - w(i, j, k - 1, ivz) * sk(i, j, k - 1, 2)
                        wwz = w(i + 1, j, k, ivz) * si(i, j, k, 3) - w(i - 1, j, k, ivz) * si(i - 1, j, k, 3) &
                              + w(i, j + 1, k, ivz) * sj(i, j, k, 3) - w(i, j - 1, k, ivz) * sj(i, j - 1, k, 3) &
                              + w(i, j, k + 1, ivz) * sk(i, j, k, 3) - w(i, j, k - 1, ivz) * sk(i, j, k - 1, 3)

                        ! (du_i/dx_j)(du_i/dx_j); the gradients above are
                        ! scaled by 2*vol.

                        fact = half / vol(i, j, k)
                        gradU2 = fact**2 * (uux**2 + uuy**2 + uuz**2 &
                                            + vvx**2 + vvy**2 + vvz**2 &
                                            + wwx**2 + wwy**2 + wwz**2)

                        ! Inverse of the local speed of sound squared.

                        a2Inv = w(i, j, k, irho) / (gamma(i, j, k) * p(i, j, k))

                        compCoef = SAc5 * a2Inv * gradU2

                        scratch(i, j, k, idvt) = scratch(i, j, k, idvt) &
                                                 - compCoef * w(i, j, k, itu1)**2

#ifndef USE_TAPENADE
                        ! -dsource/dnuTilde of this term, which is positive
                        ! for positive nuTilde.
                        qq(i, j, k) = qq(i, j, k) + max(two * compCoef * w(i, j, k, itu1), zero)
#endif
#ifdef TAPENADE_REVERSE
                    end do
#else
                end do
            end do
        end do
#endif
    end subroutine saCompressibilitySource

end module saCorrections
