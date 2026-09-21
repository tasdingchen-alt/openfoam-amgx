import sys
from paraview.simple import *

vtkFile = sys.argv[1]
rocketFile = sys.argv[2]
outPng = sys.argv[3]

reader = LegacyVTKReader(FileNames=[vtkFile])

sl = Slice(Input=reader)
sl.SliceType = 'Plane'
sl.SliceType.Origin = [0.5, 0.0, 0.0]
sl.SliceType.Normal = [0.0, 1.0, 0.0]

rocket = LegacyVTKReader(FileNames=[rocketFile])

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

# Compact scalar bar in the lower-right corner
bar = GetScalarBar(lut, view)
bar.WindowLocation = 'Any Location'
bar.Position = [0.90, 0.08]
bar.ScalarBarLength = 0.34
bar.ScalarBarThickness = 14
bar.Title = 'U (m/s)'
bar.LabelFontSize = 12
bar.TitleFontSize = 12

mdp = GetDisplayProperties(rocket, view=view)
mdp.Representation = 'Surface'
mdp.ColorArrayName = None
mdp.DiffuseColor = [0.9, 0.9, 0.9]

# Small orientation axes
view.OrientationAxesVisibility = 1
view.OrientationAxesLabelColor = [1.0, 1.0, 1.0]
view.OrientationAxesOutlineColor = [1.0, 1.0, 1.0]

Render()

view.CameraParallelProjection = 1
view.CameraFocalPoint = [0.62, 0.0, 0.0]
view.CameraPosition = [0.62, -3.0, 0.0]
view.CameraViewUp = [0.0, 0.0, 1.0]
view.CameraParallelScale = 0.6

Render()
SaveScreenshot(outPng, view)
