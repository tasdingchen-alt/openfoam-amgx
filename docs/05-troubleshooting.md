# Troubleshooting

This is a list of the problems actually hit while building and validating this
project, with their fixes.

## Build

### `fatal error: mpi.h: No such file or directory`

`PstreamGlobals.H` needs the MPI headers. `Make/options` obtains them with:

```make
$(addprefix -I,$(shell mpicc --showme:incdirs))
```

If `mpicc` is not on `PATH` when `wmake` runs, source the OpenFOAM environment
first.

### `error: 'cudaGetDeviceCount' was not declared`

`amgx_c.h` does not pull in the CUDA runtime headers. `AmgXSolver.C` includes
`<cuda_runtime.h>` explicitly.

### `cannot find -lamgx`

AmgX's shared library is **`libamgxsh.so`**, not `libamgx.so`. Link with
`-lamgxsh`.

### Blackwell / `sm_120` "unsupported gpu architecture"

The RTX 50 series needs CUDA >= 12.8. Build AmgX with
`-DCMAKE_CUDA_ARCHITECTURES=120`.

## Runtime

### Floating point exception (SIGFPE) inside AmgX

OpenFOAM v12's `etc/bashrc` contains `export FOAM_SIGFPE=` (empty but set), and
`Foam::env()` only checks whether the variable exists, so FP trapping is
enabled. AmgX performs benign FP operations that then abort the run. Always:

```bash
source /path/to/OpenFOAM/etc/bashrc
unset FOAM_SIGFPE
```

### `AMGX_solver_create` fails with `rc=1` (BAD_PARAMETERS)

The generated AmgX config is order/structure sensitive:

- `"preconditioner"` must be the **first** key of the `"solver"` object.
- The outer solver needs `"scope": "main"` and the preconditioner
  `"scope": "amg"`.
- `"monitor_residual": 0` triggers this error. Use `1`.

### `Caught amgx exception: Fail to get info from cudense`

The coarse-level dense LU found a singular (or nearly singular) coarse matrix.
For OpenFOAM's pressure matrix (row sums ~0) this happens with the default
`coarse_solver DENSE_LU_SOLVER`. Use `coarseSolver NOSOLVER` for such systems,
or make the matrix non-singular with `pRefCell`/`pRefValue`.

### Solver stagnates (final residual barely below the initial)

Almost always a wrong matrix, not an AmgX problem. Check the CSR assembly:
the halo column and coefficient conventions are subtle (see below). In this
project the bug was using `+interfaceBouCoeffs_` instead of
`-interfaceBouCoeffs_` for the processor halo entries.

### `MPI_ERR_COMM: invalid communicator`

`AMGX_resources_create` keeps a pointer to the communicator. Pass the address
of a **persistent** object, e.g. `&PstreamGlobals::MPI_COMM_FOAM`, not a local
variable.

### `D1 interpolation is not supported in distributed settings`

Classical AMG defaults to the D1 interpolator, which AmgX does not support in
parallel. Set `interpolator D2` (the solver does this automatically for
`algorithm CLASSICAL`).

### `invalid configuration argument` in `setNeighborAggregates`

A CUDA kernel was launched with zero blocks because a halo ring was empty. This
is a symptom of an incorrect halo structure (wrong global column indices),
not an AmgX bug.

### "Invalid C wrapper" / `double free or corruption` at exit

The cached AmgX handles were destroyed after AmgX's own static state was gone.
The cache deliberately does not destroy its handles at process exit.

## AmgX configuration gotchas

- Valid modes are `dDDI`, `dDFI`, `dFFI`, `hDDI`, `hDFI`, `hFFI`. There are no
  `sDDI`/`sDFI`/`sFFI` modes (single precision is `dFFI`/`hFFI`).
- `"monitor_residual"` must be `1`.
- For performance, a direct coarse solve (`DENSE_LU_SOLVER`) is often decisive
  (see `docs/04-benchmarks.md`).
