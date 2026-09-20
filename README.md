# openfoam-amgx

GPU-accelerated linear solvers for OpenFOAM, backed by NVIDIA
[AmgX](https://github.com/NVIDIA/AmgX).

`openfoam-amgx` provides a drop-in `lduMatrix::solver` that offloads OpenFOAM's
sparse linear solves to NVIDIA GPUs. It works in serial and under MPI, and is
selected from `fvSolution` like any built-in solver.

```foam
solvers
{
    p
    {
        solver      AmgX;
        AmgX        { mode dDDI; }
        tolerance   1e-6;
        relTol      0.1;
    }
}
```

## Features

- `AmgXSolver` runtime-selectable from `fvSolution`
- Serial (single GPU) and MPI parallel (multi-GPU) support
- Process-wide AmgX state cache - no rebuild on every solve
- Automatic configuration, overridable with a JSON file or dictionary options
- OpenFOAM-normalised residuals reported like the built-in solvers
- Validated against GAMG to machine precision

## How it works

OpenFOAM stores a linear system as an LDU matrix (`diag`, `upper`, `lower`) and
solves it through the `lduMatrix::solver` interface. `AmgXSolver`:

1. assembles the LDU matrix into CSR,
2. for internal face `f`: `lower[f] = A(owner, neighbour)` and
   `upper[f] = A(neighbour, owner)`,
3. in parallel, numbers cells rank-contiguously and adds the off-processor
   couplings as halo entries with global column indices,
4. hands the matrix to AmgX (`AMGX_matrix_upload_all` /
   `AMGX_matrix_upload_all_global_32`) and reads the solution back.

Boundary conditions are already folded into the diagonal and source by
OpenFOAM's `solveSegregated`, so no extra boundary handling is needed. The
pressure matrix is negative definite, so `-A`, `-b` is passed to AmgX to keep
the diagonal positive for aggregation AMG.

See `docs/00-overview.md` for details.

## Requirements

- OpenFOAM Foundation v12 (API verified against this version)
- CUDA toolkit matching your GPU (12.8+ for Blackwell `sm_120`)
- AmgX built from source, with the same MPI as OpenFOAM
- An MPI implementation matching the one OpenFOAM was built with

## Build

```bash
source /path/to/OpenFOAM/etc/bashrc
export AMGX_DIR=/path/to/amgx
export CUDA_HOME=/usr/local/cuda
./Allwmake
```

Step-by-step CUDA and AmgX installation is in `docs/02-build-amgx.md`.

## Usage

```foam
// system/controlDict
libs ("libnvidiaSolvers.so");
```

```foam
// system/fvSolution
solvers
{
    p
    {
        solver      AmgX;

        AmgX
        {
            configFile      "amgx_config.json";   // optional
            mode            dDDI;                 // optional
            smoother        MULTICOLOR_DILU;      // optional
            coarseSolver    DENSE_LU_SOLVER;      // optional (see benchmarks)
            setupEveryTime  true;                 // optional
            verbose         false;                // optional
        }

        tolerance   1e-8;
        relTol      0.01;
    }
}
```

**Important:** OpenFOAM v12 enables floating-point trapping because its
`etc/bashrc` sets `FOAM_SIGFPE=`. AmgX triggers benign FP exceptions, so always
run with:

```bash
source /path/to/OpenFOAM/etc/bashrc
unset FOAM_SIGFPE
```

Parallel:

```bash
decomposePar
mpirun -np 4 foamRun -parallel
```

## Validation and benchmarks

- Single tightly-converged `pitzDaily` pressure solve agrees with GAMG to
  machine precision (relative difference ~3e-10).
- Serial vs 2-rank parallel agree to ~2e-8 (`U`) and ~4e-9 (`p`).
- 3-D lid-driven cavity, 2 M cells: after tuning, AmgX is ~7-10% faster than
  GAMG. The **default** configuration is slower - the direct coarse solve
  (`DENSE_LU_SOLVER`) is decisive.

Details: `docs/03-validation.md` and `docs/04-benchmarks.md`.

## Documentation

| File | Contents |
|------|----------|
| `docs/00-overview.md`     | How OpenFOAM solves a system, LDU -> CSR, roadmap |
| `docs/01-quickstart.md`   | Build, usage, options, parallel |
| `docs/02-build-amgx.md`   | Installing CUDA and building AmgX |
| `docs/03-validation.md`   | Numerical validation against GAMG |
| `docs/04-benchmarks.md`   | 2 M-cell 3-D cavity benchmark vs GAMG |
| `docs/05-troubleshooting.md` | Every problem hit while building this |
| `docs/06-tutorial-stl-aero.md` | End-to-end tutorial: STL -> external aero case -> GPU solve -> animation |
| `docs/07-tutorial-cn.md` | 中文完整教程（安装 / 使用 / 实战 / 排错） |

## License

GPL-3.0-or-later (required when linking against OpenFOAM). AmgX is distributed
by NVIDIA under BSD-3-Clause.
