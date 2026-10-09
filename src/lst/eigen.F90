! Minimal LST operator export, ported verbatim from ADflowContribute src/lst/eigen.F90
! (Jacobian dR/dw and mass matrix dU/dw, PETSc binary write). Only the routines
! lst_dump.py uses; no eigensolver, preconditioner or two-level code.
module eigenapi
    use constants
    use petsc

    logical :: drdwInitialized = .False.
    logical :: massMatrixInitialized = .False.
    integer(kind=intType) :: nStateJac = 0_intType
    integer(kind=intType) :: nStateMass = 0_intType
    type(tMat) :: lstMassMat

contains

    subroutine setnstate(frozenTurbulence, nStateOut)
        use flowVarRefState, only: nwf, nw
        implicit none

        logical, intent(in) :: frozenTurbulence
        integer(kind=intType), intent(out) :: nStateOut

        ! The number of state variables differs between full and frozen
        ! turbulence formulations. We centralize this logic so Jacobian and mass
        ! matrix assembly stay consistent.
        if (frozenTurbulence) then
            nStateOut = nwf
        else
            nStateOut = nw
        end if
    end subroutine setnstate

    subroutine writematrixpetscbinary(matrixObj, fileName)
        use communication, only: ADFLOW_COMM_WORLD
        use utils, only: EChk
        implicit none

        type(tMat), intent(in) :: matrixObj
        character(len=*), intent(in) :: fileName

        type(tPetscViewer) :: viewer
        integer(kind=intType) :: ierr

        call PetscViewerBinaryOpen(ADFLOW_COMM_WORLD, trim(fileName), FILE_MODE_WRITE, viewer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call MatView(matrixObj, viewer, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        call PetscViewerDestroy(viewer, ierr)
        call EChk(ierr, __FILE__, __LINE__)
    end subroutine writematrixpetscbinary

    subroutine initializelstjacobian(usePC, frozenTurbulence)
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use flowVarRefState, only: viscous
        use stencils, only: N_visc_drdw, visc_drdw_stencil, N_euler_drdw, euler_drdw_stencil, &
                            N_visc_pc, visc_pc_stencil, N_euler_pc, euler_pc_stencil
        use ADjointPETSc, only: dRdWT
        use ADjointVars, only: nCellsLocal
        use adjointUtils, only: myMatCreate, statePreAllocation
        use utils, only: EChk
        implicit none

        logical, intent(in) :: usePC, frozenTurbulence

        integer(kind=intType) :: ierr, level, nDimW, n_stencil
        integer(kind=intType), dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
        integer(kind=intType), dimension(:, :), pointer :: stencil

        ! Re-create the Jacobian matrix if it already exists.
        if (drdwInitialized) then
            call MatDestroy(dRdWT, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            drdwInitialized = .False.
        end if

        call setnstate(frozenTurbulence, nStateJac)

        nDimW = nStateJac*nCellsLocal(1_intType)*nTimeIntervalsSpectral

        ! Preallocation follows the solver stencil because Jacobian sparsity
        ! comes from nearest-neighbor residual couplings.
        allocate (nnzDiagonal(nCellsLocal(1_intType)*nTimeIntervalsSpectral), &
                  nnzOffDiag(nCellsLocal(1_intType)*nTimeIntervalsSpectral))

        if (usePC) then
            if (viscous) then
                n_stencil = N_visc_pc
                stencil => visc_pc_stencil
            else
                n_stencil = N_euler_pc
                stencil => euler_pc_stencil
            end if
        else
            if (viscous) then
                n_stencil = N_visc_drdw
                stencil => visc_drdw_stencil
            else
                n_stencil = N_euler_drdw
                stencil => euler_drdw_stencil
            end if
        end if

        level = 1_intType
        call statePreAllocation(nnzDiagonal, nnzOffDiag, nDimW/nStateJac, stencil, n_stencil, &
                                level, .False.)
        call myMatCreate(dRdWT, nStateJac, nDimW, nDimW, nnzDiagonal, nnzOffDiag, &
                         __FILE__, __LINE__)

        ! Keep PETSc metadata consistent with expected Jacobian structure.
        call MatSetOption(dRdWT, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        deallocate (nnzDiagonal, nnzOffDiag)
        drdwInitialized = .True.
    end subroutine initializelstjacobian

    subroutine buildlstjacobian(useAD, usePC, useTranspose, useObjective, frozenTurbulence)
        use ADjointPETSc, only: dRdWT
        use adjointUtils, only: setupStateResidualMatrix
        implicit none

        logical, intent(in) :: useAD, usePC, useTranspose, useObjective, frozenTurbulence
        integer(kind=intType) :: level

        ! Lazy initialization allows a one-call build path from Python.
        if (.not. drdwInitialized) then
            call initializelstjacobian(usePC, frozenTurbulence)
        end if

        level = 1_intType
        ! Assemble the state residual Jacobian into dRdWT.
        call setupStateResidualMatrix(dRdWT, useAD, usePC, useTranspose, &
                                      useObjective, frozenTurbulence, level)
    end subroutine buildlstjacobian

    subroutine initializelstmassmatrix(frozenTurbulence)
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use ADjointVars, only: nCellsLocal
        use adjointUtils, only: myMatCreate
        use utils, only: EChk
        implicit none

        logical, intent(in) :: frozenTurbulence

        integer(kind=intType) :: ierr
        integer(kind=intType) :: nBlockRows, nDimW
        integer(kind=intType), dimension(:), allocatable :: nnzDiagonal, nnzOffDiag

        ! Re-create the mass matrix if it already exists.
        if (massMatrixInitialized) then
            call MatDestroy(lstMassMat, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            massMatrixInitialized = .False.
        end if

        call setnstate(frozenTurbulence, nStateMass)

        nBlockRows = nCellsLocal(1_intType)*nTimeIntervalsSpectral
        nDimW = nStateMass*nBlockRows

        ! dU/dw is block diagonal by cell/time-slice, so each block row only
        ! couples to itself.
        allocate (nnzDiagonal(nBlockRows), nnzOffDiag(nBlockRows))
        nnzDiagonal = 1_intType
        nnzOffDiag = 0_intType

        call myMatCreate(lstMassMat, nStateMass, nDimW, nDimW, nnzDiagonal, nnzOffDiag, &
                         __FILE__, __LINE__)

        deallocate (nnzDiagonal, nnzOffDiag)
        massMatrixInitialized = .True.
    end subroutine initializelstmassmatrix

    subroutine buildlstmassmatrix(frozenTurbulence, scale)
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use blockPointers, only: nDom, il, jl, kl, w, flowDoms
        use utils, only: EChk, setPointers, terminate
        implicit none

        logical, intent(in) :: frozenTurbulence
        real(kind=realType), intent(in) :: scale

        integer(kind=intType) :: ierr, nn, sps, i, j, k, l, iRow
#ifdef USE_COMPLEX
        complex(kind=realType), dimension(:, :), allocatable :: blk
        complex(kind=realType) :: stateScale
#else
        real(kind=realType), dimension(:, :), allocatable :: blk
        real(kind=realType) :: stateScale
#endif

        ! Lazy initialization allows direct build calls from Python.
        if (.not. massMatrixInitialized) then
            call initializelstmassmatrix(frozenTurbulence)
        end if

        if (nStateMass < 5_intType) then
            call terminate("buildlstmassmatrix", "Expected at least 5 states per cell (rho,u,v,w,rhoE).")
        end if

        ! Optional global multiplier (default 1.0). This mirrors the previous
        ! Python mass-matrix utility behavior.
        stateScale = scale

        allocate (blk(nStateMass, nStateMass))
        ! Preserve sparsity structure/allocation, but overwrite all values.
        call MatZeroEntries(lstMassMat, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1_intType, sps)
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            iRow = flowDoms(nn, 1_intType, sps)%globalCell(i, j, k)

                            blk = 0.0_realType
                            ! For non-compute cells, keep an identity block so
                            ! the matrix remains well-defined at those entries.
                            if (flowDoms(nn, 1_intType, sps)%iBlank(i, j, k) .ne. 1_intType) then
                                do l = 1, nStateMass
                                    blk(l, l) = stateScale
                                end do
                            else
                                ! Conservative-vs-state map:
                                ! U = [rho, rho*u, rho*v, rho*w, rhoE, turb...]
                                ! w = [rho, u, v, w, rhoE, turb...]
                                ! so dU/dw has rho-velocity cross terms and rho
                                ! on momentum diagonals.
                                blk(1, 1) = stateScale
                                blk(2, 1) = stateScale*w(i, j, k, 2)
                                blk(2, 2) = stateScale*w(i, j, k, 1)
                                blk(3, 1) = stateScale*w(i, j, k, 3)
                                blk(3, 3) = stateScale*w(i, j, k, 1)
                                blk(4, 1) = stateScale*w(i, j, k, 4)
                                blk(4, 4) = stateScale*w(i, j, k, 1)
                                blk(5, 5) = stateScale
                                if (nStateMass >= 6_intType) then
                                    do l = 6, nStateMass
                                        blk(l, l) = stateScale
                                    end do
                                end if
                            end if

                            call MatSetValuesBlocked(lstMassMat, 1, iRow, 1, iRow, blk, INSERT_VALUES, ierr)
                            call EChk(ierr, __FILE__, __LINE__)
                        end do
                    end do
                end do
            end do
        end do

        ! Finalize PETSc assembly so downstream operations (save/use) are valid.
        call MatAssemblyBegin(lstMassMat, MAT_FINAL_ASSEMBLY, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        call MatAssemblyEnd(lstMassMat, MAT_FINAL_ASSEMBLY, ierr)
        call EChk(ierr, __FILE__, __LINE__)

        deallocate (blk)
    end subroutine buildlstmassmatrix

    subroutine writelstjacobian(fileName, negate)
        use adjointAPI, only: saveADjointMatrix
        use ADjointPETSc, only: dRdWT
        use utils, only: EChk, terminate
        implicit none

        character(len=*), intent(in) :: fileName
        logical, intent(in) :: negate
        integer(kind=intType) :: ierr
#ifdef USE_COMPLEX
        complex(kind=realType) :: negOne
#else
        real(kind=realType) :: negOne
#endif

        ! Guard against writing an uninitialized matrix handle.
        if (.not. drdwInitialized) then
            call terminate("writelstjacobian", "drdw has not been initialized.")
        end if

        negOne = -1.0_realType

        ! Optional sign flip is done in-place before write, then restored.
        ! This avoids duplicating matrix storage just for sign convention.
        if (negate) then
            call MatScale(dRdWT, negOne, ierr)
            call EChk(ierr, __FILE__, __LINE__)
        end if

        call saveADjointMatrix(trim(fileName))

        if (negate) then
            call MatScale(dRdWT, negOne, ierr)
            call EChk(ierr, __FILE__, __LINE__)
        end if
    end subroutine writelstjacobian

    subroutine writelstmassmatrix(fileName)
        use utils, only: terminate
        implicit none

        character(len=*), intent(in) :: fileName

        ! Guard against writing an uninitialized matrix handle.
        if (.not. massMatrixInitialized) then
            call terminate("writelstmassmatrix", "dudw has not been initialized.")
        end if

        call writematrixpetscbinary(lstMassMat, trim(fileName))
    end subroutine writelstmassmatrix

    subroutine destroylstjacobian
        use ADjointPETSc, only: dRdWT
        use utils, only: EChk
        implicit none

        integer(kind=intType) :: ierr

        ! Release PETSc storage after Jacobian export path.
        if (drdwInitialized) then
            call MatDestroy(dRdWT, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            drdwInitialized = .False.
        end if
        nStateJac = 0_intType
    end subroutine destroylstjacobian

    subroutine destroylstmassmatrix
        use utils, only: EChk
        implicit none

        integer(kind=intType) :: ierr

        ! Release PETSc storage after mass-matrix export path.
        if (massMatrixInitialized) then
            call MatDestroy(lstMassMat, ierr)
            call EChk(ierr, __FILE__, __LINE__)
            massMatrixInitialized = .False.
        end if
        nStateMass = 0_intType
    end subroutine destroylstmassmatrix

end module eigenapi
