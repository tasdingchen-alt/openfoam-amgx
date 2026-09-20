#!/bin/bash
# Build and install AmgX into $INSTALL_PREFIX (default: $HOME/amgx).
#
# Usage:
#   CUDA_HOME=/usr/local/cuda-12-8 ./install-amgx.sh        # sm_120 (Blackwell)
#   CUDA_HOME=/usr/local/cuda-12-8 ./install-amgx.sh 90     # Hopper

set -e

CUDA_HOME="${CUDA_HOME:-/usr/local/cuda}"
ARCH="${1:-120}"
INSTALL_PREFIX="${INSTALL_PREFIX:-$HOME/amgx}"
SRC_DIR="${SRC_DIR:-$HOME/AmgX}"
JOBS="$(nproc)"

sudo apt-get install -y build-essential cmake git

if [ ! -d "$SRC_DIR" ]; then
    git clone https://github.com/NVIDIA/AmgX.git "$SRC_DIR"
fi

cd "$SRC_DIR"
rm -rf build
mkdir build
cd build

cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" \
    -DCMAKE_CUDA_ARCHITECTURES="$ARCH" \
    -DCMAKE_CUDA_COMPILER="$CUDA_HOME/bin/nvcc"

make -j"$JOBS"
make install

echo
echo "AmgX installed to $INSTALL_PREFIX"
ls -l "$INSTALL_PREFIX/lib/libamgxsh.so" "$INSTALL_PREFIX/include/amgx_c.h"
