#!/bin/bash
# Install the CUDA toolkit inside WSL (Ubuntu).
#
# Do NOT install a driver inside WSL - the Windows driver is shared. Only the
# toolkit is needed. Use the "wsl-ubuntu" repository for this reason.
#
# Usage:
#   ./install-cuda-wsl.sh                 # install CUDA 12.8 toolkit
#   ./install-cuda-wsl.sh cuda-toolkit-13-0

set -e

CUDA_PKG="${1:-cuda-toolkit-12-8}"

sudo apt-get update
sudo apt-get install -y wget gnupg

KEYRING=/tmp/cuda-keyring_1.1-1_all.deb
wget -O "$KEYRING" \
    https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i "$KEYRING"
sudo apt-get update

sudo apt-get install -y "$CUDA_PKG"

echo
echo "Installed. Add the following to ~/.bashrc and re-source it:"
echo
echo "    export CUDA_HOME=/usr/local/cuda-12-8"
echo "    export PATH=\$CUDA_HOME/bin:\$PATH"
echo "    export LD_LIBRARY_PATH=\$CUDA_HOME/lib64:\$LD_LIBRARY_PATH"
echo
echo "Then check with: nvcc --version"
