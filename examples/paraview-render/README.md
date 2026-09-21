# ParaView 渲染参数与镜头设置（AIM-120D / 水火箭）

本文档记录两个气动算例的**动画渲染镜头参数**，以及本项目在 ParaView 上踩过的坑。
以后渲染这两个模型时，直接照抄即可，不要随意改镜头。

---

## 一、必须知道的 4 个 ParaView 坑

1. **`Render()` 会自动重置相机**
   第一次 `Render()` 会把相机自动适配到整个数据（`CameraParallelScale` 被改掉）。
   → **必须"先 Render 一次，再设相机，再 Render"**，否则相机设置无效。

2. **`CameraParallelScale` 语义与直觉相反、且非线性**
   数值越**小**画面越**近**。本项目实测（1280×720，正交投影）：
   | scale | 效果 |
   |---|---|
   | 0.04 | 场铺满整个画面 |
   | 0.20 | 场铺满（偏近） |
   | 0.60 | 场约占画面 52% 宽（**水火箭采用**） |
   | 1.2  | 视野很大（AIM-120D 采用） |
   → 调镜头时**只改这一个数**，按"再近/再远"微调，不要动其它相机参数。

3. **色标窗口位置必须写 `'Any Location'`**（带空格），写成 `'AnyLocation'` 会报
   `obsolete value` 错误。

4. **动画每帧必须用独立 `pvbatch` 进程渲染**
   `foamToVTK` 写的 VTK 无 TIME 字段，ParaView 会把时间步当 0..N 且渲染有缓存，
   同一进程内切帧会得到完全相同的图。
   → 用 shell 循环，**每帧启动一个 `pvbatch`**，读单个 VTK 文件渲染。

---

## 二、AIM-120D（沿 z 轴，机头在 -z）

- 几何：`AIM120D_scaled.stl`，长 3.66 m，**机头在 z=-1.83，尾翼在 z=+1.83**
- 来流：沿 +z（从机头吹向机尾），入口在 z 小端

### 渲染参数（侧视，机头朝左）
```python
# 切片：对称面 y=0
sl.SliceType = 'Plane'
sl.SliceType.Origin = [0.0, 0.0, 0.3]
sl.SliceType.Normal = [0.0, 1.0, 0.0]

# 相机（在第一次 Render() 之后设置）
view.CameraParallelProjection = 1
view.CameraFocalPoint = [0.0, 0.0, 0.3]
view.CameraPosition   = [0.0, 9.0, 0.3]
view.CameraViewUp     = [1.0, 0.0, 0.0]     # x 朝上 → 导弹长轴(z)水平
view.CameraParallelScale = 1.2              # ← 镜头远近，只调这个
```
- 物理量：速度大小 `U`（Jet，0.33–120）或压力 `p`（Cool to Warm，−10000~6000）
- 分辨率：1920×800 或 1280×720

---

## 三、水火箭（沿 x 轴，机头在 x=0）

- 几何：`water-rocket-straightened-smooth.stl`，长 1.243 m，直径 ~0.154 m
  - **机头在 x=0**（半径小），**尾翼在 x≈1.2**（半径大）
  - 曲面非封闭、法向不一致（snappyHexMesh 仍可正常贴合）
- 来流：沿 +x（从机头吹向机尾），入口在 x 小端，30 m/s

### 渲染参数（侧视，机头朝左）——**已定稿，保持不变**
```python
# 切片：对称面 y=0
sl.SliceType = 'Plane'
sl.SliceType.Origin = [0.5, 0.0, 0.0]
sl.SliceType.Normal = [0.0, 1.0, 0.0]

# 相机（在第一次 Render() 之后设置）
view.CameraParallelProjection = 1
view.CameraFocalPoint = [0.62, 0.0, 0.0]
view.CameraPosition   = [0.62, -3.0, 0.0]
view.CameraViewUp     = [0.0, 0.0, 1.0]     # z 朝上 → 火箭长轴(x)水平
view.CameraParallelScale = 0.6              # ★ 定稿值，不要改
```
- 物理量：速度大小 `U`（Jet，**0.33–120**）
- 分辨率：**1280×720（720p）**
- 色标：紧凑竖条，右下角
  ```python
  bar = GetScalarBar(lut, view)
  bar.WindowLocation = 'Any Location'
  bar.Position = [0.90, 0.08]
  bar.ScalarBarLength = 0.34
  bar.ScalarBarThickness = 14
  bar.Title = 'U (m/s)'
  bar.LabelFontSize = 12
  bar.TitleFontSize = 12
  ```
- 坐标轴：白色小号，左下角
  ```python
  view.OrientationAxesVisibility = 1
  view.OrientationAxesLabelColor = [1.0, 1.0, 1.0]
  view.OrientationAxesOutlineColor = [1.0, 1.0, 1.0]
  ```

---

## 四、渲染脚本模板

### 单帧脚本 `render_rocket.py`（参数：场文件、表面文件、输出 png）
```python
import sys
from paraview.simple import *

vtkFile = sys.argv[1]
surfaceFile = sys.argv[2]
outPng = sys.argv[3]

reader = LegacyVTKReader(FileNames=[vtkFile])

sl = Slice(Input=reader)
sl.SliceType = 'Plane'
sl.SliceType.Origin = [0.5, 0.0, 0.0]
sl.SliceType.Normal = [0.0, 1.0, 0.0]

surface = LegacyVTKReader(FileNames=[surfaceFile])

view = GetActiveViewOrCreate('RenderView')
view.ViewSize = [1280, 720]
view.BackgroundColorMode = 'Single Color'
view.Background = [0.05, 0.05, 0.08]

dp = GetDisplayProperties(sl, view=view)
ColorBy(dp, ('POINTS', 'U', 'Magnitude'))
lut = GetColorTransferFunction('U')
lut.RescaleTransferFunction(0.33, 120.0)
lut.ApplyPreset('Jet', True)
dp.SetScalarBarVisibility(view, True)
bar = GetScalarBar(lut, view)
bar.WindowLocation = 'Any Location'
bar.Position = [0.90, 0.08]
bar.ScalarBarLength = 0.34
bar.ScalarBarThickness = 14
bar.Title = 'U (m/s)'

mdp = GetDisplayProperties(surface, view=view)
mdp.Representation = 'Surface'
mdp.ColorArrayName = None
mdp.DiffuseColor = [0.9, 0.9, 0.9]

# --- 关键：先 Render 一次，再设相机 ---
Render()
view.CameraParallelProjection = 1
view.CameraFocalPoint = [0.62, 0.0, 0.0]
view.CameraPosition = [0.62, -3.0, 0.0]
view.CameraViewUp = [0.0, 0.0, 1.0]
view.CameraParallelScale = 0.6
Render()
SaveScreenshot(outPng, view)
```

### 逐帧脚本 `renderall.sh`（每帧独立进程）
```bash
#!/bin/bash
cd /tmp/rocket || exit 1
mkdir -p frames
rm -f frames/*.png
i=0
for f in $(ls -1 VTK/rocket_*.vtk | sort -V); do
    step=$(basename "$f" | sed 's/rocket_//;s/\.vtk//')
    pvbatch render_rocket.py "$f" "VTK/rocket/rocket_${step}.vtk" \
        "$(printf 'frames/frame_%04d.png' $i)" > /dev/null 2>&1
    i=$((i + 1))
done
echo RENDERED
ls frames | wc -l
```

### 合成视频
```bash
ffmpeg -y -framerate 12 -i frames/frame_%04d.png \
    -c:v libx264 -crf 18 -pix_fmt yuv420p out.mp4
ffmpeg -y -framerate 12 -i frames/frame_%04d.png \
    -vf 'fps=12,scale=960:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse' \
    out.gif
```

---

## 五、参数速查

| 项 | AIM-120D | 水火箭 |
|---|---|---|
| 长轴 | z | x |
| 机头位置 | z = −1.83 | x = 0 |
| 流向 | +z | +x |
| 切片原点 | [0, 0, 0.3] | [0.5, 0, 0] |
| 切片法向 | [0, 1, 0] | [0, 1, 0] |
| 相机位置 | [0, 9, 0.3] | [0.62, −3, 0] |
| 相机焦点 | [0, 0, 0.3] | [0.62, 0, 0] |
| 相机上向 | [1, 0, 0] | [0, 0, 1] |
| **CameraParallelScale** | **1.2** | **0.6（定稿）** |
| 分辨率 | 1920×800 / 1280×720 | **1280×720** |
| 物理量 | U 或 p | U（0.33–120, Jet） |
