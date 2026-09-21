#!/bin/bash
cd /tmp/rocket || exit 1
mkdir -p frames
rm -f frames/*.png
i=0
for f in $(ls -1 VTK/rocket_*.vtk 2>/dev/null | sort -V); do
    step=$(basename "$f" | sed 's/rocket_//;s/\.vtk//')
    rf="VTK/rocket/rocket_${step}.vtk"
    pvbatch render_rocket.py "$f" "$rf" "$(printf 'frames/frame_%04d.png' $i)" > /dev/null 2>&1
    i=$((i + 1))
done
echo RENDERED
ls frames | wc -l
