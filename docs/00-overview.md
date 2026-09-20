# Overview

## Goal

Provide GPU-accelerated linear solvers for OpenFOAM by implementing a new
`lduMatrix::solver` that assembles OpenFOAM's LDU matrix into a CSR matrix and
hands it to NVIDIA AmgX.

## How OpenFOAM solves a linear system

1. `fvMatrix<Type>::solve()` calls `solveSegregated()` (or `solveCoupled()`).
2. `solveSegregated()`:
   - copies `source_` and calls `addBoundarySource(source)`,
   - for each component: copies `diag()`, applies `addBoundaryDiag(diag(), cmpt)`,
     extracts `boundaryCoeffs_` / `internalCoeffs_` and the interface list,
   - corrects the source for the explicit part of coupled boundaries via
     `initMatrixInterfaces` / `updateMatrixInterfaces`,
   - constructs a solver with `lduMatrix::solver::New(...)` and calls
     `solve(psiCmpt, sourceCmpt, cmpt)`.
3. The solver receives an `lduMatrix` view (`diag`, `upper`, `lower`) plus the
   interface coefficients and interface list.

## LDU to CSR

For internal face `f` with `owner = upperAddr[f]` and `neighbour = lowerAddr[f]`:

- `lower[f]` is `A(owner, neighbour)` and goes into row `owner`, column `neighbour`.
- `upper[f]` is `A(neighbour, owner)` and goes into row `neighbour`, column `owner`.
- `diag[i]` is `A(i, i)` and already contains the implicit boundary contribution.

Coupled boundaries (processor, cyclic) contribute additional off-diagonal
entries whose column is the neighbour cell index obtained from the interface
addressing. Non-coupled boundaries contribute only to `diag` and `source`.

## Status

| Phase | Deliverable                                   | Status |
|-------|-----------------------------------------------|--------|
| 0     | Toolchain and AmgX build                      | done   |
| 1     | Serial scalar `AmgXSolver`                    | done   |
| 2     | Process-wide AmgX state cache                 | done   |
| 3     | MPI / multi-GPU (halo columns)                | done   |
| 4     | Dictionary options and tolerance mapping      | done   |
| 5     | CPU fallback on failure                       | todo   |
| 6     | Validation and benchmarks                     | done   |
| 7     | Packaging, documentation, GitHub release      | done   |

Validation: `docs/03-validation.md`. Benchmarks: `docs/04-benchmarks.md`.
