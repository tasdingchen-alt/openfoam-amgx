# Tutorial: external aerodynamics from an STL with AmgX

This is an end-to-end walkthrough: take a missile/vehicle STL, build an
external-aerodynamics case in OpenFOAM Foundation v12, solve the pressure
equation with `openfoam-amgx` on the GPU, and render an animation.

It is written so that a person **or an AI assistant** can follow it
step by step. Every file and command is given in full.

The example geometry is `AIM120D_scaled.stl` (a slender body along the z-axis),
but the recipe works for any closed STL.

## 0. Prerequisites

- OpenFOAM Foundation v12 (`source /path/to/OpenFOAM/etc/bashrc`)
- CUDA + AmgX built as in `docs/02-build-amgx.md`
- `openfoam-amgx` built (see the main README)
- ParaView (`pvbatch`) and `ffmpeg` for the animation

**Always** `unset FOAM_SIGFPE` before running (see `docs/05-troubleshooting.md`).

## 1. Inspect the STL and find the orientation

The most common mistake is getting the flow direction wrong. First find where
the nose is. Save this as `stlprof.py` and run it:

```python
import struct
import numpy as np

path = 'constant/triSurface/body.stl'
with open(path, 'rb') as f:
    f.read(80)
    n = struct.unpack('<I', f.read(4))[0]
    tris = np.zeros((n, 3, 3), dtype=float)
    for i in range(n):
        d = struct.unpack('<12fH', f.read(50))
        tris[i, 0], tris[i, 1], tris[i, 2] = d[3:6], d[6:9], d[9:12]

v = tris.reshape(-1, 3)
z = v[:, 2]
r = np.sqrt(v[:, 0]**2 + v[:, 1]**2)

nb = 20
edges = np.linspace(z.min(), z.max(), nb + 1)
for i in range(nb):
    m = (z >= edges[i]) & (z < edges[i + 1])
    if m.sum():
        print('%7.3f  r_max=%7.4f' % (0.5*(edges[i]+edges[i+1]), r[m].max()))
```

The **nose is the end with the small radius**; the **tail (fins) has the large
radius**. In our case the nose is at `z ≈ -1.83` and the fins at `z ≈ +1.83`.

**Rule:** the freestream must hit the nose. If the nose is at `-z`, the inlet
goes on the `-z` face with `U = (0, 0, +50)`.

## 2. Case layout

```
case/
  0/U 0/p 0/k 0/omega 0/nut
  constant/physicalProperties
  constant/momentumTransport
  constant/triSurface/body.stl
  system/blockMeshDict
  system/snappyHexMeshDict
  system/controlDict  fvSchemes  fvSolution
```

Create it and copy the STL:

```bash
mkdir -p case/0 case/constant/triSurface case/system
cp AIM120D_scaled.stl case/constant/triSurface/body.stl
cd case
```

## 3. Background mesh

`system/blockMeshDict` — a box around the body, longer downstream. Nose at `-z`,
so the inlet is at `z = -5` and the outlet at `z = +8`:

```
FoamFile { version 2.0; format ascii; class dictionary; object blockMeshDict; }

scale 1;

vertices
(
    (-2 -2 -5) ( 2 -2 -5) ( 2  2 -5) (-2  2 -5)
    (-2 -2  8) ( 2 -2  8) ( 2  2  8) (-2  2  8)
);

blocks ( hex (0 1 2 3 4 5 6 7) (20 20 65) simpleGrading (1 1 1) );

boundary
(
    inlet  { type patch; faces ((0 3 2 1)); }   // z = -5  (upstream, nose side)
    outlet { type patch; faces ((4 5 6 7)); }   // z = +8  (downstream)
    walls  { type patch; faces ((0 4 7 3)(1 5 4 0)(2 6 5 1)(3 7 6 2)); }
);

mergePatchPairs ();
```

```bash
blockMesh
```

## 4. Snap the mesh to the STL

`system/snappyHexMeshDict`. Note the **v12 geometry syntax** uses `file`:

```
castellatedMesh true;
snap            true;
addLayers       true;

geometry
{
    missile
    {
        type triSurfaceMesh;
        file "body.stl";
    }
}

castellatedMeshControls
{
    maxLocalCells       2000000;
    maxGlobalCells      4000000;
    minRefinementCells  10;
    nCellsBetweenLevels 3;

    features ();

    refinementSurfaces
    {
        missile { level (3 4); patchInfo { type wall; } }
    }

    resolveFeatureAngle 30;
    refinementRegions   {}
    locationInMesh (1.5 1.5 0);
    allowFreeStandingZoneFaces true;
}

snapControls
{
    nSmoothPatch 3; tolerance 2.0; nSolveIter 100;
    nRelaxIter 5; nFeatureSnapIter 10;
}

addLayersControls
{
    relativeSizes true;
    layers { missile { nSurfaceLayers 3; } }
    expansionRatio 1.2; finalLayerThickness 0.5; minThickness 0.1;
    nGrow 0; featureAngle 60; slipFeatureAngle 30;
    nRelaxIter 3; nSmoothSurfaceNormals 1; nSmoothNormals 3;
    nSmoothThickness 10; maxFaceThicknessRatio 0.5;
    maxThicknessToMedialRatio 0.3; minMedianAxisAngle 90;
    nBufferCellsNoExtrude 0; nLayerIter 50;
}

meshQualityControls
{
    maxNonOrtho 65; maxBoundarySkewness 20; maxInternalSkewness 4;
    maxConcave 80; minVol 1e-13; minTetQuality 1e-15;
    minArea -1; minTwist 0.02; minDeterminant 0.001;
    minFaceWeight 0.02; minVolRatio 0.01; minTriangleTwist -1;
    nSmoothScale 4; errorReduction 0.75;
}

debug 0;
mergeTolerance 1e-6;
```

```bash
snappyHexMesh -overwrite
checkMesh | grep -E 'cells:|Mesh OK'
```

Expected: a few hundred thousand cells and `Mesh OK`. The body appears as the
`missile` patch.

## 5. Fields

`0/U` (flow along +z, hitting the nose at `-z`):

```
FoamFile { format ascii; class volVectorField; object U; }
dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 50);
boundaryField
{
    inlet   { type fixedValue;   value uniform (0 0 50); }
    outlet  { type zeroGradient; }
    walls   { type slip; }
    missile { type noSlip; }
}
```

`0/p`:

```
FoamFile { format ascii; class volScalarField; object p; }
dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
    inlet   { type zeroGradient; }
    outlet  { type fixedValue; value uniform 0; }
    walls   { type zeroGradient; }
    missile { type zeroGradient; }
}
```

`0/k` (`k = 1.5 (I U)^2`, `I = 0.05`), `0/omega` (`~1000`), `0/nut`:

```
// k
internalField uniform 9.375;
boundaryField
{
    inlet   { type fixedValue; value uniform 9.375; }
    outlet  { type zeroGradient; }
    walls   { type zeroGradient; }
    missile { type kqRWallFunction; value uniform 9.375; }
}

// omega
internalField uniform 1000;
boundaryField
{
    inlet   { type fixedValue; value uniform 1000; }
    outlet  { type zeroGradient; }
    walls   { type zeroGradient; }
    missile { type omegaWallFunction; value uniform 1000; }
}

// nut
internalField uniform 0;
boundaryField
{
    inlet   { type calculated; value uniform 0; }
    outlet  { type calculated; value uniform 0; }
    walls   { type calculated; value uniform 0; }
    missile { type nutkWallFunction; value uniform 0; }
}
```

`constant/physicalProperties`:

```
viscosityModel constant;
nu 1.5e-05;
```

`constant/momentumTransport`:

```
simulationType RAS;
RAS { model kOmegaSST; turbulence on; printCoeffs on; }
```

## 6. Solver settings

`system/fvSchemes`:

```
ddtSchemes { default backward; }
gradSchemes { default Gauss linear; grad(U) cellLimited Gauss linear 1; }
divSchemes
{
    default none;
    div(phi,U)     Gauss linearUpwind grad(U);
    div(phi,k)     Gauss upwind;
    div(phi,omega) Gauss upwind;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
wallDist { method meshWave; }
```

`system/fvSolution` — **this is where AmgX is selected**:

```
solvers
{
    p
    {
        solver      AmgX;
        tolerance   1e-7;
        relTol      0.01;

        AmgX
        {
            mode            dDDI;
            smoother        MULTICOLOR_DILU;
            coarseSolver    DENSE_LU_SOLVER;   // decisive for performance
            maxLevels       25;
        }
    }
    pFinal { $p; relTol 0; }

    "(U|k|omega)"
    {
        solver      smoothSolver;
        smoother    symGaussSeidel;
        tolerance   1e-7;
        relTol      0.1;
    }
    "(U|k|omega)Final" { $U; relTol 0; }
}

PIMPLE
{
    nOuterCorrectors 1;
    nCorrectors      2;
    nNonOrthogonalCorrectors 1;
}
```

`system/controlDict`:

```
application     foamRun;
solver          incompressibleFluid;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         0.07;
deltaT          0.0005;
writeControl    adjustableRunTime;
writeInterval   0.002;
writeFormat     binary;
runTimeModifiable no;
adjustTimeStep  yes;
maxCo           5;
libs            ("libnvidiaSolvers.so");
```

## 7. Run

```bash
source /path/to/OpenFOAM/etc/bashrc
unset FOAM_SIGFPE                 # essential
export AMGX_DIR=/path/to/amgx
export CUDA_HOME=/usr/local/cuda
foamRun | tee log.run
```

Look for lines like:

```
AmgX:  Solving for p, Initial residual = 1, Final residual = ..., No Iterations 9
```

## 8. Animation

Convert to VTK and render one PNG per time step with ParaView, then encode.

`foamToVTK` writes `VTK/body_<step>.vtk` plus `VTK/<patch>/<patch>_<step>.vtk`.

`render1.py` (renders one frame; arguments: field file, body file, output png):

```python
import sys
from paraview.simple import *

reader = LegacyVTKReader(FileNames=[sys.argv[1]])
sl = Slice(Input=reader)
sl.SliceType = 'Plane'
sl.SliceType.Origin = [0.0, 0.0, 0.3]
sl.SliceType.Normal = [0.0, 1.0, 0.0]      # symmetry plane y = 0

body = LegacyVTKReader(FileNames=[sys.argv[2]])

view = GetActiveViewOrCreate('RenderView')
view.ViewSize = [1280, 720]
view.BackgroundColorMode = 'Single Color'
view.Background = [0.05, 0.05, 0.08]
view.CameraParallelProjection = 1
view.CameraFocalPoint = [0.0, 0.0, 0.3]
view.CameraPosition = [0.0, 9.0, 0.3]
view.CameraViewUp = [1.0, 0.0, 0.0]
view.CameraParallelScale = 1.2

dp = GetDisplayProperties(sl, view=view)
ColorBy(dp, ('POINTS', 'p'))
lut = GetColorTransferFunction('p')
lut.RescaleTransferFunction(-900.0, 700.0)
lut.ApplyPreset('Cool to Warm', True)
dp.SetScalarBarVisibility(view, True)

mdp = GetDisplayProperties(body, view=view)
mdp.Representation = 'Surface'
mdp.ColorArrayName = None
mdp.DiffuseColor = [0.9, 0.9, 0.9]

Render()
SaveScreenshot(sys.argv[3], view)
```

Render every frame in a separate process (avoids ParaView time-step caching):

```bash
#!/bin/bash
foamToVTK
mkdir -p frames
i=0
for f in $(ls -1 VTK/*_*.vtk | grep -v /missile/ | sort -V); do
    step=$(basename "$f" | sed 's/.*_//;s/\.vtk//')
    pvbatch render1.py "$f" "VTK/missile/missile_${step}.vtk" \
        "$(printf 'frames/frame_%04d.png' $i)" > /dev/null 2>&1
    i=$((i + 1))
done
ffmpeg -y -framerate 12 -i frames/frame_%04d.png \
    -c:v libx264 -crf 18 -pix_fmt yuv420p animation.mp4
```

## 9. Checklist

- [ ] STL orientation checked (nose upstream)
- [ ] `unset FOAM_SIGFPE`
- [ ] `checkMesh` reports `Mesh OK`
- [ ] `libs ("libnvidiaSolvers.so")` in `controlDict`
- [ ] `solver AmgX` for `p` in `fvSolution`
- [ ] `DENSE_LU_SOLVER` coarse solver for speed (see `docs/04-benchmarks.md`)

## 10. Common pitfalls

| Symptom | Fix |
|---------|-----|
| SIGFPE inside AmgX | `unset FOAM_SIGFPE` |
| Wake appears upstream of the body | nose/tail swapped - see section 1 |
| `keyword file is undefined` in `snappyHexMeshDict` | v12 needs `type triSurfaceMesh; file "x.stl";` |
| `surfaceFeatureExtract` says it is superseded | use `surfaceFeatures` in v12, or drop `features` |
| Animation frames all identical | render each frame in a separate `pvbatch` process |
| AmgX slower than GAMG | set `coarseSolver DENSE_LU_SOLVER`; see `docs/04-benchmarks.md` |

See `docs/05-troubleshooting.md` for the full list.
