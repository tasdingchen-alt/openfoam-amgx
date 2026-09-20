# Installing CUDA and AmgX

Verified against AmgX `main` (2025), which uses `CMAKE_CUDA_ARCHITECTURES` and
already lists `120` (Blackwell) among the allowed architectures.

Target environment for this guide: **WSL2 Ubuntu 24.04**, RTX 50-series
(Blackwell, `sm_120`), OpenFOAM Foundation v12 using system OpenMPI 4.1.6.

## 1. Check the GPU

```bash
nvidia-smi
```

If this works inside WSL, the Windows driver is being passed through
correctly. Do **not** install a driver inside WSL.

## 2. Install the CUDA toolkit (WSL)

Use the `wsl-ubuntu` repository, which ships the toolkit without a driver.
CUDA **12.8 or newer** is required for `sm_120`.

```bash
sudo apt-get update
sudo apt-get install -y wget gnupg

wget https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
sudo apt-get update

sudo apt-get install -y cuda-toolkit-12-8
```

Add to `~/.bashrc`:

```bash
export CUDA_HOME=/usr/local/cuda-12-8
export PATH=$CUDA_HOME/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH
```

Then:

```bash
source ~/.bashrc
nvcc --version
```

The scripts `scripts/install-cuda-wsl.sh` automates the above.

## 3. Build AmgX

```bash
sudo apt-get install -y build-essential cmake git

git clone https://github.com/NVIDIA/AmgX.git ~/AmgX
cd ~/AmgX
mkdir -p build && cd build

cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$HOME/amgx \
    -DCMAKE_CUDA_ARCHITECTURES=120 \
    -DCMAKE_CUDA_COMPILER=$CUDA_HOME/bin/nvcc

make -j$(nproc)
make install
```

`-DCMAKE_CUDA_ARCHITECTURES=120` narrows the build to Blackwell only, which
is much faster. Omit it to build the default `90;100;120`.

The script `scripts/install-amgx.sh` automates the above.

Verify:

```bash
ls $HOME/amgx/include/amgx_c.h
ls $HOME/amgx/lib/libamgxsh.so      # note: amgxsh, not amgx
```

## 4. Build openfoam-amgx

Keep the source on the Linux filesystem (not `/mnt/c`) for build speed:

```bash
cp -r "/mnt/c/Users/<you>/.../openfoam-amgx" ~/openfoam-amgx

source /opt/openfoam12/etc/bashrc
export AMGX_DIR=$HOME/amgx
export CUDA_HOME=/usr/local/cuda-12-8
export LD_LIBRARY_PATH=$AMGX_DIR/lib:$CUDA_HOME/lib64:$LD_LIBRARY_PATH

cd ~/openfoam-amgx
./Allwmake
```

The library is installed to `$FOAM_USER_LIBBIN/libnvidiaSolvers.so`.

## Notes

- **MPI**: OpenFOAM v12 here uses `openmpi-system` (system OpenMPI 4.1.6).
  AmgX's `find_package(MPI)` picks up the same system OpenMPI, so they match.
  A mismatch causes hangs in parallel runs.
- **Library name**: the AmgX shared library is `libamgxsh.so`, hence
  `-lamgxsh` in `Make/options`.
- **CUDA 13.1+**: AmgX switches to C++17 automatically.
- **No GPU visible**: ensure `/usr/lib/wsl/lib` is on `LD_LIBRARY_PATH`.
