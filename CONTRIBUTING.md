# Contributing

Thanks for your interest in `openfoam-amgx`!

## Reporting issues

Please include:

- OpenFOAM version (`foamVersion` / `foamSystemCheck` output)
- CUDA version (`nvcc --version`) and GPU
- AmgX version / commit
- The relevant `fvSolution` / `controlDict` entries
- The full error output

Many problems are already covered in `docs/05-troubleshooting.md` - please check
there first.

## Development

The library is built with `wmake`:

```bash
source /path/to/OpenFOAM/etc/bashrc
export AMGX_DIR=/path/to/amgx
export CUDA_HOME=/usr/local/cuda
./Allwmake
```

The code targets OpenFOAM Foundation v12. The public API used is
`lduMatrix::solver`; porting to other versions mainly means adjusting that
interface.

## Pull requests

- Keep changes focused and match the existing style.
- Update the relevant file under `docs/` if behaviour changes.
- Describe how you tested the change (case, mesh size, serial/parallel).

## Code layout

```
src/nvidiaSolvers/
  Make/                      wmake files and options
  AmgXSolver/
    AmgXSolver.H/.C          the lduMatrix::solver implementation
    ldu2csr.H/.C             LDU -> CSR assembly (with halo support)
```

## License

By contributing you agree that your contributions are licensed under
GPL-3.0-or-later, consistent with the project.
