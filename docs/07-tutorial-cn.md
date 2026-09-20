# openfoam-amgx 中文使用教程

本教程带你从零开始：安装环境 → 编译 `openfoam-amgx` → 在 `fvSolution` 里使用
GPU 求解器 → 从 STL 几何做外部气动算例 → 出动画。每一步都有可直接复制的命令
和完整文件内容，也适合作为 AI 助手的配置参考。

---

## 一、这是什么

OpenFOAM 的线性方程组（压力、动量等）默认用 CPU 求解器（GAMG、PCG…）。
`openfoam-amgx` 提供一个可运行时选择的 `lduMatrix::solver`，把线性求解交给
NVIDIA GPU（通过 [AmgX](https://github.com/NVIDIA/AmgX) 库）。

用法和内置求解器完全一样，只改 `fvSolution`：

```foam
p
{
    solver      AmgX;
    AmgX        { mode dDDI; }
    tolerance   1e-7;
    relTol      0.01;
}
```

**原理**：OpenFOAM 的矩阵是 LDU 格式（`diag/upper/lower`）。求解器把它组装成
CSR，并行时按 rank 连续编号、把跨进程耦合作为 halo 列（全局列号）传给 AmgX，
解完再写回。边界条件已被 OpenFOAM 折进对角和源项，无需额外处理。

---

## 二、环境准备

目标环境：Ubuntu（含 WSL2）、OpenFOAM Foundation v12、NVIDIA GPU。

### 1. 检查 GPU

```bash
nvidia-smi
```

WSL 里能看到 GPU 即可。**不要在 WSL 里装显卡驱动**，驱动由 Windows 透传。

### 2. 安装 CUDA Toolkit（12.8 及以上）

RTX 50 系（Blackwell, sm_120）**必须 CUDA ≥ 12.8**。用 NVIDIA 的 WSL 源：

```bash
sudo apt-get update
sudo apt-get install -y wget gnupg

wget https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
sudo apt-get update
sudo apt-get install -y cuda-toolkit-12-8
```

写入 `~/.bashrc`：

```bash
export CUDA_HOME=/usr/local/cuda-12.8
export PATH=$CUDA_HOME/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH
```

```bash
source ~/.bashrc
nvcc --version
```

### 3. 编译 AmgX

```bash
sudo apt-get install -y build-essential cmake git

git clone https://github.com/NVIDIA/AmgX.git ~/AmgX
cd ~/AmgX && mkdir -p build && cd build

cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$HOME/amgx \
    -DCMAKE_CUDA_ARCHITECTURES=120 \
    -DCMAKE_CUDA_COMPILER=$CUDA_HOME/bin/nvcc

make -j$(nproc)
make install
```

验证（注意是 **`libamgxsh.so`**）：

```bash
ls $HOME/amgx/include/amgx_c.h $HOME/amgx/lib/libamgxsh.so
```

> AmgX 用系统 MPI 编译（`find_package(MPI)`），与 OpenFOAM 用的 MPI 一致即可。

### 4. 确认 OpenFOAM

```bash
ls /opt/openfoam12/etc/bashrc
source /opt/openfoam12/etc/bashrc
foamVersion
```

---

## 三、编译 openfoam-amgx

```bash
# 建议放在 Linux 文件系统里（不要放 /mnt/c，编译慢）
cp -r openfoam-amgx ~/openfoam-amgx

source /opt/openfoam12/etc/bashrc
export AMGX_DIR=$HOME/amgx
export CUDA_HOME=/usr/local/cuda-12.8

cd ~/openfoam-amgx
./Allwmake
```

成功后得到 `$FOAM_USER_LIBBIN/libnvidiaSolvers.so`。

---

## 四、使用（串行）

### 1. 加载库

`system/controlDict`：

```foam
libs ("libnvidiaSolvers.so");
```

### 2. 在 fvSolution 里选 AmgX

```foam
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
            coarseSolver    DENSE_LU_SOLVER;   // 性能关键，见基准文档
            maxLevels       25;
        }
    }
    pFinal { $p; relTol 0; }
}
```

### 3. 运行前必须做的一件事

OpenFOAM v12 的 `etc/bashrc` 里有 `export FOAM_SIGFPE=`（空但已定义），会误开
浮点异常捕获，AmgX 的正常浮点运算会触发 SIGFPE 崩溃。**运行前必须**：

```bash
source /opt/openfoam12/etc/bashrc
unset FOAM_SIGFPE
foamRun
```

---

## 五、使用（并行 MPI）

```bash
unset FOAM_SIGFPE
decomposePar
mpirun -np 4 foamRun -parallel
```

- 每个 rank 分配一个 GPU：rank `r` 用 `r % 显卡数`。单卡机器上所有 rank 共用
  device 0（结果正确，但小算例比串行慢）。
- AmgX 与 OpenFOAM 必须用**同一套 MPI**。
- 已验证：串行与 2 rank 并行结果一致（U 相对差 2e-8，p 4e-9）。

---

## 六、配置选项一览（`AmgX` 子字典）

| 选项 | 默认 | 说明 |
|---|---|---|
| `configFile` | 无 | AmgX JSON 配置文件，给定后覆盖下列所有项 |
| `mode` | `dDDI` | 精度/索引模式：`dDDI dDFI dFFI hDDI hDFI hFFI` |
| `verbose` | `false` | 打印 AmgX 求解统计 |
| `setupEveryTime` | `true` | 每次求解都重建 AMG 层级（安全） |
| `solver` | `FGMRES` | `FGMRES GMRES PCG PBICGSTAB AMG` |
| `algorithm` | `AGGREGATION` | `AGGREGATION CLASSICAL` |
| `smoother` | `BLOCK_JACOBI` | 如 `MULTICOLOR_DILU`、`MULTICOLOR_GS` |
| `coarseSolver` | `NOSOLVER` | `NOSOLVER` 或 `DENSE_LU_SOLVER` |
| `interpolator` | `D2` | 仅 CLASSICAL 用；D1 并行不支持 |
| `selector` | `SIZE_2` | 聚合选择器 |
| `maxLevels` | `50` | AMG 最大层数 |

**容差映射**：`relTol>0` 时用相对容差 `relTol`；否则用绝对 `tolerance`。
OpenFOAM 归一化后的初始残差约为 1，所以两者接近。

**性能提示**：`coarseSolver DENSE_LU_SOLVER`（直接粗解）对性能影响最大，
迭代数可降数倍。但若矩阵奇异（纯 Neumann 压力方程），需改用 `NOSOLVER`，
或用 `pRefCell/pRefValue` 设参考单元。

---

## 七、完整实战：从 STL 到动画

以导弹 STL 为例（细长体沿 z 轴）。

### 1. 判断机头方向（最关键，最容易错）

用下面脚本看沿 z 的半径分布，**半径小的一端是机头，大的一端是尾翼**：

```python
# stlprof.py
import struct, numpy as np
path = 'constant/triSurface/body.stl'
with open(path, 'rb') as f:
    f.read(80)
    n = struct.unpack('<I', f.read(4))[0]
    tris = np.zeros((n, 3, 3))
    for i in range(n):
        d = struct.unpack('<12fH', f.read(50))
        tris[i,0], tris[i,1], tris[i,2] = d[3:6], d[6:9], d[9:12]
v = tris.reshape(-1, 3); z = v[:,2]
r = np.sqrt(v[:,0]**2 + v[:,1]**2)
edges = np.linspace(z.min(), z.max(), 21)
for i in range(20):
    m = (z >= edges[i]) & (z < edges[i+1])
    if m.sum():
        print('%7.3f  r_max=%7.4f' % (0.5*(edges[i]+edges[i+1]), r[m].max()))
```

**规则：来流必须吹向机头。** 若机头在 `-z`，则入口在 `-z` 面，来流
`U = (0,0,+50)`。

### 2. 目录结构

```
case/
  0/U 0/p 0/k 0/omega 0/nut
  constant/physicalProperties  constant/momentumTransport
  constant/triSurface/body.stl
  system/blockMeshDict  snappyHexMeshDict
  system/controlDict  fvSchemes  fvSolution
```

### 3. 背景网格 `system/blockMeshDict`

机头在 `-z` → 入口在 `z=-5`，出口在 `z=+8`：

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
    inlet  { type patch; faces ((0 3 2 1)); }   // z=-5，机头侧（上游）
    outlet { type patch; faces ((4 5 6 7)); }   // z=+8（下游）
    walls  { type patch; faces ((0 4 7 3)(1 5 4 0)(2 6 5 1)(3 7 6 2)); }
);

mergePatchPairs ();
```

```bash
blockMesh
```

### 4. 贴合网格 `system/snappyHexMeshDict`

注意 v12 的 geometry 用 **`file`** 关键字：

```
castellatedMesh true;
snap            true;
addLayers       true;

geometry
{
    missile { type triSurfaceMesh; file "body.stl"; }
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

### 5. 初始/边界场

`0/U`（来流 +z，打机头）：

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

`0/p`：

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

`0/k` / `0/omega` / `0/nut`：

```
// k：k = 1.5 (I U)^2，I=0.05，U=50 → 9.375
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

`constant/physicalProperties`：

```
viscosityModel constant;
nu 1.5e-05;
```

`constant/momentumTransport`：

```
simulationType RAS;
RAS { model kOmegaSST; turbulence on; printCoeffs on; }
```

### 6. 求解器设置

`system/fvSchemes`：

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

`system/fvSolution`（AmgX 在这里启用）：

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
            coarseSolver    DENSE_LU_SOLVER;
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

`system/controlDict`：

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

### 7. 运行

```bash
source /opt/openfoam12/etc/bashrc
unset FOAM_SIGFPE
export AMGX_DIR=$HOME/amgx
export CUDA_HOME=/usr/local/cuda-12.8
foamRun | tee log.run
```

日志里会出现：

```
AmgX:  Solving for p, Initial residual = 1, Final residual = ..., No Iterations 9
```

### 8. 出动画

```bash
foamToVTK
```

`foamToVTK` 会写 `VTK/body_<步>.vtk` 和 `VTK/<patch>/<patch>_<步>.vtk`。

`render1.py`（单帧渲染，参数：场文件、弹体文件、输出 png）：

```python
import sys
from paraview.simple import *

reader = LegacyVTKReader(FileNames=[sys.argv[1]])
sl = Slice(Input=reader)
sl.SliceType = 'Plane'
sl.SliceType.Origin = [0.0, 0.0, 0.3]
sl.SliceType.Normal = [0.0, 1.0, 0.0]      # 对称面 y=0

body = LegacyVTKReader(FileNames=[sys.argv[2]])

view = GetActiveViewOrCreate('RenderView')
view.ViewSize = [1280, 720]
view.BackgroundColorMode = 'Single Color'
view.Background = [0.05, 0.05, 0.08]
view.CameraParallelProjection = 1
view.CameraFocalPoint = [0.0, 0.0, 0.3]
view.CameraPosition = [0.0, 9.0, 0.3]
view.CameraViewUp = [1.0, 0.0, 0.0]
view.CameraParallelScale = 1.2             # 数字越小，画面越放大

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

逐帧渲染脚本（**每帧一个独立进程**，避免 ParaView 时间步缓存）：

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

# 合成 mp4 和 gif
ffmpeg -y -framerate 12 -i frames/frame_%04d.png \
    -c:v libx264 -crf 18 -pix_fmt yuv420p animation.mp4
ffmpeg -y -framerate 12 -i frames/frame_%04d.png \
    -vf 'fps=12,scale=720:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse' \
    animation.gif
```

---

## 八、结果验证与性能

- **数值正确性**：紧容差单步，AmgX 与 GAMG 结果一致到机器精度
  （p 相对差 ~3e-10）；串行与并行一致（U 2e-8）。
- **性能**（3D 方腔，200 万单元，10 步）：

  | 配置 | 时间 |
  |---|---|
  | GAMG | 96 s |
  | AmgX 默认（`NOSOLVER` 粗解） | 154 s（更慢） |
  | AmgX 调优（`MULTICOLOR_DILU` + `DENSE_LU_SOLVER`） | **87 s** |

  **结论**：默认配置可能比 GAMG 慢，**`DENSE_LU_SOLVER` 直接粗解是关键**。
  GPU 不一定更快，务必在自己的算例上实测调参。

---

## 九、常见问题排查

| 现象 | 解决 |
|---|---|
| AmgX 内部触发 SIGFPE | 运行前 `unset FOAM_SIGFPE` |
| 尾流出现在机头前方 | 机头/尾翼方向反了，见第七节第 1 步 |
| `AMGX_solver_create` 报 rc=1 | 配置需 `preconditioner` 在最前、`scope: main/amg`、`monitor_residual: 1` |
| `Fail to get info from cudense` | 粗网格矩阵奇异，改用 `coarseSolver NOSOLVER` 或设 `pRefCell` |
| 求解停滞（残差几乎不降） | halo 列/系数错了，检查并行组装 |
| `MPI_ERR_COMM: invalid communicator` | 传给 AmgX 的 communicator 必须是持久对象 |
| `D1 interpolation is not supported` | CLASSICAL 要设 `interpolator D2` |
| 退出时 `double free` | AmgX 句柄不能在进程退出时才销毁（本项目已处理） |
| `snappyHexMeshDict` 报 `keyword file is undefined` | v12 用 `type triSurfaceMesh; file "x.stl";` |
| `surfaceFeatureExtract` 提示已被取代 | v12 改用 `surfaceFeatures`，或去掉 `features` |
| 动画每帧都一样 | 每帧用独立 `pvbatch` 进程渲染 |
| 链接找不到 `-lamgx` | 库名是 `libamgxsh.so`，用 `-lamgxsh` |
| Blackwell 报不支持架构 | 需要 CUDA ≥ 12.8，AmgX 用 `-DCMAKE_CUDA_ARCHITECTURES=120` |

更完整的排查见 `docs/05-troubleshooting.md`。

---

## 十、文件清单

| 文件 | 内容 |
|---|---|
| `docs/00-overview.md` | 原理与状态 |
| `docs/01-quickstart.md` | 快速开始、选项表、并行 |
| `docs/02-build-amgx.md` | CUDA/AmgX 安装 |
| `docs/03-validation.md` | 数值验证 |
| `docs/04-benchmarks.md` | 性能基准 |
| `docs/05-troubleshooting.md` | 全部踩坑 |
| `docs/06-tutorial-stl-aero.md` | 英文版 STL 气动教程 |
| `docs/07-tutorial-cn.md` | 本中文教程 |

许可证：GPL-3.0-or-later（链接 OpenFOAM 的要求）；AmgX 为 BSD-3-Clause。
