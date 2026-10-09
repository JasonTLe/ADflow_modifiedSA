! This module contains the source and diffusion terms of variants of the
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
!   3. The density correction of the diffusion term of Catris and Aupoix
!      (SA-Catris), selected with densityCorrectionSA = "Catris-Aupoix"
!      (original, non-conservative form) or "Catris-Aupoix-cons"
!      (conservative form), see saCatrisViscous and
!      saDensityCorrectionSource:
!
!        Catris, S. and Aupoix, B., "Density Corrections for Turbulence
!        Models," Aerospace Science and Technology, Vol. 4, No. 1, 2000,
!        pp. 1-11, https://doi.org/10.1016/S1270-9638(00)00112-7.
!
!   4. The density-gradient diffusion terms of the conservative form of
!      the SA equation ("conservative"), of Catris (1999, "Catris-1999")
!      and of Seror et al. (J. Aircraft 42(1), 2005, "Seror"), see
!      saDensityCorrectionSource and Xue et al., AIAA J. 62(1), 2024,
!      https://doi.org/10.2514/1.J063274, for a summary.
!
! All corrections can be combined with each other and with "SA" or
! "SA-Edwards".
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
        !       correction if useCompressibilitySA is set and the density
        !       correction selected by densityCorrSA.
        !
        use constants
        use blockPointers, only: il, jl, kl
        use inputPhysics, only: turbModel, useCompressibilitySA, densityCorrSA
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

        if (densityCorrSA /= saDensityNone .and. densityCorrSA /= saDensityCatrisAupoixCons) then
            call saDensityCorrectionSource
        end if

        ! Advection Term
        nn = itu1 - 1
        call turbAdvection(1_intType, 1_intType, nn, qq)

        ! Unsteady Term
        call unsteadyTurbTerm(1_intType, 1_intType, nn, qq)

        ! Viscous Terms
        if (densityCorrSA == saDensityCatrisAupoix .or. densityCorrSA == saDensityCatrisAupoixCons) then
            call saCatrisViscous
        else
            call saViscous
        end if

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

    subroutine saCatrisViscous
        !
        !  Diffusion term of the Catris-Aupoix density-corrected SA model,
        !  which replaces saViscous when densityCorrSA == saDensityCatrisAupoix:
        !
        !    (1/sigma) [ d/dx_j (mu dnuTilde/dx_j)
        !                + d/dx_j (phi dphi/dx_j) + cb2 (dphi/dx_j)^2 ] / rho,
        !    phi = sqrt(rho) nuTilde,
        !
        !  Catris, S. and Aupoix, B., "Density Corrections for Turbulence
        !  Models," Aerospace Science and Technology, Vol. 4, No. 1, 2000,
        !  pp. 1-11, https://doi.org/10.1016/S1270-9638(00)00112-7
        !  (SA-Catris on the NASA Turbulence Modeling Resource). The equation
        !  for rho nuTilde is divided by rho to obtain the non-conservative
        !  form solved by ADflow. The discretization follows saViscous:
        !  cross derivatives are neglected and the terms in phi are written
        !  as div((1+cb2) phi grad(phi)) - cb2 phi div(grad(phi)) with
        !  nonnegative coefficients. For constant density the term reduces
        !  to the standard SA diffusion term.

        use blockPointers
        use constants
        use paramTurb
#ifndef USE_TAPENADE
        use sa, only: qq
#endif
        implicit none

        ! Local variables.
        integer(kind=intType) :: i, j, k, ii
        real(kind=realType) :: cb3Inv, rhoInv
        real(kind=realType) :: voli, volmi, volpi, xm, ym, zm, xp, yp, zp
        real(kind=realType) :: xa, ya, za, ttm, ttp
        real(kind=realType) :: sqr0, sqrm, sqrp, phi0, phim, phip
        real(kind=realType) :: clm, clp, ctm, ctp
        real(kind=realType) :: c1m, c1p, c10, b1, c1, d1

        cb3Inv = one / rsaCb3

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

                        ! phi = sqrt(rho) nuTilde in the cell and at the faces.

                        sqr0 = sqrt(w(i, j, k, irho))
                        sqrm = sqrt(w(i, j, k - 1, irho))
                        sqrp = sqrt(w(i, j, k + 1, irho))
                        phi0 = sqr0 * w(i, j, k, itu1)
                        phim = half * (sqrm * w(i, j, k - 1, itu1) + phi0)
                        phip = half * (sqrp * w(i, j, k + 1, itu1) + phi0)

                        ! Laminar diffusion coefficients, mu/sigma at the faces.

                        clm = half * (rlv(i, j, k - 1) + rlv(i, j, k)) * ttm * cb3Inv
                        clp = half * (rlv(i, j, k + 1) + rlv(i, j, k)) * ttp * cb3Inv

                        ! Turbulent diffusion coefficients multiplying the
                        ! differences of phi, ((1+cb2) phi_face - cb2 phi)/sigma,
                        ! see saViscous. As there, the total coefficients of the
                        ! neighbouring values (laminar + turbulent) are kept
                        ! nonnegative, which is needed at walls where the halo
                        ! value of nuTilde is negative.

                        ctm = max(((one + rsaCb2) * phim - rsaCb2 * phi0) * ttm * cb3Inv, -clm / sqrm)
                        ctp = max(((one + rsaCb2) * phip - rsaCb2 * phi0) * ttp * cb3Inv, -clp / sqrp)

                        ! Coefficients of nuTilde in the neighbours and the cell
                        ! itself; the whole term is divided by rho.

                        rhoInv = one / w(i, j, k, irho)
                        c1m = (clm + ctm * sqrm) * rhoInv
                        c1p = (clp + ctp * sqrp) * rhoInv
                        c10 = (clm + clp + (ctm + ctp) * sqr0) * rhoInv

                        scratch(i, j, k, idvt) = scratch(i, j, k, idvt) + c1m * w(i, j, k - 1, itu1) &
                                                 - c10 * w(i, j, k, itu1) + c1p * w(i, j, k + 1, itu1)
#ifndef USE_TAPENADE
                        ! Update the central jacobian, see saViscous.

                        b1 = -c1m
                        c1 = c10
                        d1 = -c1p

                        if (k == 2) then
                            qq(i, j, k) = qq(i, j, k) + c1 &
                                          - b1 * max(bmtk1(i, j, itu1, itu1), zero)
                        else if (k == kl) then
                            qq(i, j, k) = qq(i, j, k) + c1 &
                                          - d1 * max(bmtk2(i, j, itu1, itu1), zero)
                        else
                            qq(i, j, k) = qq(i, j, k) + c1
                        end if
#endif
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

                        ! phi = sqrt(rho) nuTilde in the cell and at the faces.

                        sqr0 = sqrt(w(i, j, k, irho))
                        sqrm = sqrt(w(i, j - 1, k, irho))
                        sqrp = sqrt(w(i, j + 1, k, irho))
                        phi0 = sqr0 * w(i, j, k, itu1)
                        phim = half * (sqrm * w(i, j - 1, k, itu1) + phi0)
                        phip = half * (sqrp * w(i, j + 1, k, itu1) + phi0)

                        ! Laminar diffusion coefficients, mu/sigma at the faces.

                        clm = half * (rlv(i, j - 1, k) + rlv(i, j, k)) * ttm * cb3Inv
                        clp = half * (rlv(i, j + 1, k) + rlv(i, j, k)) * ttp * cb3Inv

                        ! Turbulent diffusion coefficients multiplying the
                        ! differences of phi, ((1+cb2) phi_face - cb2 phi)/sigma,
                        ! see saViscous. As there, the total coefficients of the
                        ! neighbouring values (laminar + turbulent) are kept
                        ! nonnegative, which is needed at walls where the halo
                        ! value of nuTilde is negative.

                        ctm = max(((one + rsaCb2) * phim - rsaCb2 * phi0) * ttm * cb3Inv, -clm / sqrm)
                        ctp = max(((one + rsaCb2) * phip - rsaCb2 * phi0) * ttp * cb3Inv, -clp / sqrp)

                        ! Coefficients of nuTilde in the neighbours and the cell
                        ! itself; the whole term is divided by rho.

                        rhoInv = one / w(i, j, k, irho)
                        c1m = (clm + ctm * sqrm) * rhoInv
                        c1p = (clp + ctp * sqrp) * rhoInv
                        c10 = (clm + clp + (ctm + ctp) * sqr0) * rhoInv

                        scratch(i, j, k, idvt) = scratch(i, j, k, idvt) + c1m * w(i, j - 1, k, itu1) &
                                                 - c10 * w(i, j, k, itu1) + c1p * w(i, j + 1, k, itu1)
#ifndef USE_TAPENADE
                        ! Update the central jacobian, see saViscous.

                        b1 = -c1m
                        c1 = c10
                        d1 = -c1p

                        if (j == 2) then
                            qq(i, j, k) = qq(i, j, k) + c1 &
                                          - b1 * max(bmtj1(i, k, itu1, itu1), zero)
                        else if (j == jl) then
                            qq(i, j, k) = qq(i, j, k) + c1 &
                                          - d1 * max(bmtj2(i, k, itu1, itu1), zero)
                        else
                            qq(i, j, k) = qq(i, j, k) + c1
                        end if
#endif
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

                        ! phi = sqrt(rho) nuTilde in the cell and at the faces.

                        sqr0 = sqrt(w(i, j, k, irho))
                        sqrm = sqrt(w(i - 1, j, k, irho))
                        sqrp = sqrt(w(i + 1, j, k, irho))
                        phi0 = sqr0 * w(i, j, k, itu1)
                        phim = half * (sqrm * w(i - 1, j, k, itu1) + phi0)
                        phip = half * (sqrp * w(i + 1, j, k, itu1) + phi0)

                        ! Laminar diffusion coefficients, mu/sigma at the faces.

                        clm = half * (rlv(i - 1, j, k) + rlv(i, j, k)) * ttm * cb3Inv
                        clp = half * (rlv(i + 1, j, k) + rlv(i, j, k)) * ttp * cb3Inv

                        ! Turbulent diffusion coefficients multiplying the
                        ! differences of phi, ((1+cb2) phi_face - cb2 phi)/sigma,
                        ! see saViscous. As there, the total coefficients of the
                        ! neighbouring values (laminar + turbulent) are kept
                        ! nonnegative, which is needed at walls where the halo
                        ! value of nuTilde is negative.

                        ctm = max(((one + rsaCb2) * phim - rsaCb2 * phi0) * ttm * cb3Inv, -clm / sqrm)
                        ctp = max(((one + rsaCb2) * phip - rsaCb2 * phi0) * ttp * cb3Inv, -clp / sqrp)

                        ! Coefficients of nuTilde in the neighbours and the cell
                        ! itself; the whole term is divided by rho.

                        rhoInv = one / w(i, j, k, irho)
                        c1m = (clm + ctm * sqrm) * rhoInv
                        c1p = (clp + ctp * sqrp) * rhoInv
                        c10 = (clm + clp + (ctm + ctp) * sqr0) * rhoInv

                        scratch(i, j, k, idvt) = scratch(i, j, k, idvt) + c1m * w(i - 1, j, k, itu1) &
                                                 - c10 * w(i, j, k, itu1) + c1p * w(i + 1, j, k, itu1)
#ifndef USE_TAPENADE
                        ! Update the central jacobian, see saViscous.

                        b1 = -c1m
                        c1 = c10
                        d1 = -c1p

                        if (i == 2) then
                            qq(i, j, k) = qq(i, j, k) + c1 &
                                          - b1 * max(bmti1(j, k, itu1, itu1), zero)
                        else if (i == il) then
                            qq(i, j, k) = qq(i, j, k) + c1 &
                                          - d1 * max(bmti2(j, k, itu1, itu1), zero)
                        else
                            qq(i, j, k) = qq(i, j, k) + c1
                        end if
#endif
#ifdef TAPENADE_REVERSE
                    end do
#else
                end do
            end do
        end do
#endif
    end subroutine saCatrisViscous

    subroutine saDensityCorrectionSource
        !
        !  Source terms of the density corrections of the SA model that are
        !  not part of saCatrisViscous, written in the non-conservative
        !  form for nuTilde solved by ADflow.
        !
        !  densityCorrSA == saDensityCatrisAupoix: the original Catris-Aupoix
        !  model has D(rho nuTilde)/Dt on the left hand side (TMR:
        !  SA-noft2-Catris), which differs from rho D(nuTilde)/Dt by
        !  nuTilde D(rho)/Dt = -rho nuTilde div(u). The term nuTilde div(u)
        !  is therefore added; without it (saDensityCatrisAupoixCons, TMR:
        !  SA-noft2-CatrisCons, Waligura et al. AIAA 2022-0587) the model
        !  is conservative in rho nuTilde.
        !
        !  densityCorrSA == saDensityConservative, saDensityCatris1999 or
        !  saDensitySeror: density-gradient diffusion terms. ADflow's
        !  standard diffusion term, multiplied by rho, equals the
        !  conservative form
        !
        !    (1/sigma) [ d/dx_j (rho (nu + nuTilde) dnuTilde/dx_j)
        !                + cb2 rho (dnuTilde/dx_j)^2 ]
        !
        !  minus the cross term (1/sigma) (nu + nuTilde) (drho/dx_j)(dnuTilde/dx_j)
        !  (Allmaras, Johnson and Spalart, ICCFD7-1902, 2012). The following
        !  terms are therefore added to the source, divided by rho:
        !
        !    conservative:  (1/sigma) (nu + nuTilde) grad(rho).grad(nuTilde)
        !    Catris-1999:   conservative + (cb2/sigma) nuTilde grad(rho).grad(nuTilde)
        !                   (Catris, PhD thesis, ENSAE 1999; Xue et al., AIAA J.
        !                   62(1), 2024, Eq. 5)
        !    Seror:         Catris-1999 + (cb2/sigma) grad(mu).grad(nuTilde)
        !                   (Seror et al., J. Aircraft 42(1), 2005)
        !
        !  All terms are treated explicitly (no contribution to qq).

        use blockPointers
        use constants
        use paramTurb
        use inputPhysics, only: densityCorrSA
        implicit none

        ! Local variables.
        integer(kind=intType) :: i, j, k, ii
        real(kind=realType) :: fact, rhox, rhoy, rhoz, nux, nuy, nuz, mux, muy, muz
        real(kind=realType) :: gradRhoNu, nu, coef, src, cb3Inv, divU

        cb3Inv = one / rsaCb3

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
                        fact = half / vol(i, j, k)

                        if (densityCorrSA == saDensityCatrisAupoix) then

                            ! Original (non-conservative) form of the
                            ! Catris-Aupoix model, D(rho nuTilde)/Dt on the left
                            ! hand side: divided by rho this adds nuTilde div(u)
                            ! to the equation for nuTilde solved here.

                            divU = fact * (w(i + 1, j, k, ivx) * si(i, j, k, 1) - w(i - 1, j, k, ivx) * si(i - 1, j, k, 1) &
                                           + w(i, j + 1, k, ivx) * sj(i, j, k, 1) - w(i, j - 1, k, ivx) * sj(i, j - 1, k, 1) &
                                           + w(i, j, k + 1, ivx) * sk(i, j, k, 1) - w(i, j, k - 1, ivx) * sk(i, j, k - 1, 1) &
                                           + w(i + 1, j, k, ivy) * si(i, j, k, 2) - w(i - 1, j, k, ivy) * si(i - 1, j, k, 2) &
                                           + w(i, j + 1, k, ivy) * sj(i, j, k, 2) - w(i, j - 1, k, ivy) * sj(i, j - 1, k, 2) &
                                           + w(i, j, k + 1, ivy) * sk(i, j, k, 2) - w(i, j, k - 1, ivy) * sk(i, j, k - 1, 2) &
                                           + w(i + 1, j, k, ivz) * si(i, j, k, 3) - w(i - 1, j, k, ivz) * si(i - 1, j, k, 3) &
                                           + w(i, j + 1, k, ivz) * sj(i, j, k, 3) - w(i, j - 1, k, ivz) * sj(i, j - 1, k, 3) &
                                           + w(i, j, k + 1, ivz) * sk(i, j, k, 3) - w(i, j, k - 1, ivz) * sk(i, j, k - 1, 3))

                            src = w(i, j, k, itu1) * divU

                        else

                            ! Cell-centred gradients of rho, nuTilde and mu
                            ! (Green-Gauss with the neighbouring cell values,
                            ! see saSource).

                            rhox = fact * (w(i + 1, j, k, irho) * si(i, j, k, 1) - w(i - 1, j, k, irho) * si(i - 1, j, k, 1) &
                                           + w(i, j + 1, k, irho) * sj(i, j, k, 1) - w(i, j - 1, k, irho) * sj(i, j - 1, k, 1) &
                                           + w(i, j, k + 1, irho) * sk(i, j, k, 1) - w(i, j, k - 1, irho) * sk(i, j, k - 1, 1))
                            rhoy = fact * (w(i + 1, j, k, irho) * si(i, j, k, 2) - w(i - 1, j, k, irho) * si(i - 1, j, k, 2) &
                                           + w(i, j + 1, k, irho) * sj(i, j, k, 2) - w(i, j - 1, k, irho) * sj(i, j - 1, k, 2) &
                                           + w(i, j, k + 1, irho) * sk(i, j, k, 2) - w(i, j, k - 1, irho) * sk(i, j, k - 1, 2))
                            rhoz = fact * (w(i + 1, j, k, irho) * si(i, j, k, 3) - w(i - 1, j, k, irho) * si(i - 1, j, k, 3) &
                                           + w(i, j + 1, k, irho) * sj(i, j, k, 3) - w(i, j - 1, k, irho) * sj(i, j - 1, k, 3) &
                                           + w(i, j, k + 1, irho) * sk(i, j, k, 3) - w(i, j, k - 1, irho) * sk(i, j, k - 1, 3))

                            nux = fact * (w(i + 1, j, k, itu1) * si(i, j, k, 1) - w(i - 1, j, k, itu1) * si(i - 1, j, k, 1) &
                                          + w(i, j + 1, k, itu1) * sj(i, j, k, 1) - w(i, j - 1, k, itu1) * sj(i, j - 1, k, 1) &
                                          + w(i, j, k + 1, itu1) * sk(i, j, k, 1) - w(i, j, k - 1, itu1) * sk(i, j, k - 1, 1))
                            nuy = fact * (w(i + 1, j, k, itu1) * si(i, j, k, 2) - w(i - 1, j, k, itu1) * si(i - 1, j, k, 2) &
                                          + w(i, j + 1, k, itu1) * sj(i, j, k, 2) - w(i, j - 1, k, itu1) * sj(i, j - 1, k, 2) &
                                          + w(i, j, k + 1, itu1) * sk(i, j, k, 2) - w(i, j, k - 1, itu1) * sk(i, j, k - 1, 2))
                            nuz = fact * (w(i + 1, j, k, itu1) * si(i, j, k, 3) - w(i - 1, j, k, itu1) * si(i - 1, j, k, 3) &
                                          + w(i, j + 1, k, itu1) * sj(i, j, k, 3) - w(i, j - 1, k, itu1) * sj(i, j - 1, k, 3) &
                                          + w(i, j, k + 1, itu1) * sk(i, j, k, 3) - w(i, j, k - 1, itu1) * sk(i, j, k - 1, 3))

                            gradRhoNu = rhox * nux + rhoy * nuy + rhoz * nuz

                            ! (nu + nuTilde) grad(rho).grad(nuTilde): conservative form

                            nu = rlv(i, j, k) / w(i, j, k, irho)
                            coef = nu + w(i, j, k, itu1)

                            ! + cb2 nuTilde grad(rho).grad(nuTilde): Catris (1999)

                            if (densityCorrSA == saDensityCatris1999 .or. densityCorrSA == saDensitySeror) then
                                coef = coef + rsaCb2 * w(i, j, k, itu1)
                            end if

                            src = coef * gradRhoNu

                            ! + cb2 grad(mu).grad(nuTilde): Seror et al. (2005)

                            if (densityCorrSA == saDensitySeror) then
                                mux = fact * (rlv(i + 1, j, k) * si(i, j, k, 1) - rlv(i - 1, j, k) * si(i - 1, j, k, 1) &
                                              + rlv(i, j + 1, k) * sj(i, j, k, 1) - rlv(i, j - 1, k) * sj(i, j - 1, k, 1) &
                                              + rlv(i, j, k + 1) * sk(i, j, k, 1) - rlv(i, j, k - 1) * sk(i, j, k - 1, 1))
                                muy = fact * (rlv(i + 1, j, k) * si(i, j, k, 2) - rlv(i - 1, j, k) * si(i - 1, j, k, 2) &
                                              + rlv(i, j + 1, k) * sj(i, j, k, 2) - rlv(i, j - 1, k) * sj(i, j - 1, k, 2) &
                                              + rlv(i, j, k + 1) * sk(i, j, k, 2) - rlv(i, j, k - 1) * sk(i, j, k - 1, 2))
                                muz = fact * (rlv(i + 1, j, k) * si(i, j, k, 3) - rlv(i - 1, j, k) * si(i - 1, j, k, 3) &
                                              + rlv(i, j + 1, k) * sj(i, j, k, 3) - rlv(i, j - 1, k) * sj(i, j - 1, k, 3) &
                                              + rlv(i, j, k + 1) * sk(i, j, k, 3) - rlv(i, j, k - 1) * sk(i, j, k - 1, 3))
                                src = src + rsaCb2 * (mux * nux + muy * nuy + muz * nuz)
                            end if

                            ! The terms are divided by sigma and rho.

                            src = cb3Inv * src / w(i, j, k, irho)

                        end if

                        scratch(i, j, k, idvt) = scratch(i, j, k, idvt) + src
#ifdef TAPENADE_REVERSE
                    end do
#else
                end do
            end do
        end do
#endif
    end subroutine saDensityCorrectionSource

end module saCorrections
