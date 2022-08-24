fastSHT is a very fast toolkit for doing spherical harmonic transforms on a large number of spherical maps. It converts massive SHT operations to a BLAS level 3 problem and uses the highly optimized matrix multiplication toolkit to accelerate the computation. GPU acceleration is supported and can be very effective. Core code is written in fortran but a Python wrapper is provided and recommended.


# Dependencies

Fortran compiler: `ifort` is recommanded for the CPU version; `nvfortran` is required for the GPU version

Intel MKL library

[`f90wrap`](https://github.com/jameskermode/f90wrap)

`Python3`, `numpy`, `CMake`

# Installation


## Before compilation 

Dependencies and system PATH should be configured properly. 

### (Recommended) Download and compile with the docker image compatable with both CPU and GPU version

See https://docs.docker.com/engine/install/ for a docker installation instruction

With docker installed, use

```
sudo docker pull rectaflex/intel_nvidia_sdk
```
to pull the docker image.

To enable GPU in a docker container, the Nvidia container runtime is needed, and it can be installed by 

```
curl -s -L https://nvidia.github.io/nvidia-container-runtime/gpgkey |   sudo apt-key add -
distribution=$(. /etc/os-release;echo $ID$VERSION_ID)
curl -s -L https://nvidia.github.io/nvidia-container-runtime/$distribution/nvidia-container-runtime.list |   sudo tee /etc/apt/sources.list.d/nvidia-container-runtime.list
sudo apt-get update
sudo apt-get install nvidia-container-runtime
sudo systemctl restart docker
```
Finally, run the docker image with
```
sudo docker run -it -v /home:/home --gpus all rectaflex/intel_nvidia_sdk
```

### Environment preparation for ubuntu

For non-docker users, we give a sample script for building a compilation environment for fastSHT based on an ubuntu-20.04.

```
sudo apt update

sudo apt -y install pkg-config sudo build-essential

# installing recent cmake
sudo apt-get install apt-transport-https ca-certificates gnupg software-properties-common wget

wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | sudo apt-key add -

sudo apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main'
sudo apt-get update

sudo apt-get install cmake

# Installing intel-one-api

wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
| gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null

echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list

sudo apt-get update

sudo apt install intel-basekit
sudo apt install intel-hpckit

sed -i '1 i\source /opt/intel/oneapi/setvars.sh > /dev/null' ~/.bashrc

# installing nvidia hpc sdk

echo 'deb [trusted=yes] https://developer.download.nvidia.com/hpc-sdk/ubuntu/amd64 /' | sudo tee /etc/apt/sources.list.d/nvhpc.list

sudo apt-get update -y

sudo apt-get install -y nvhpc-22-7

sed -i '1 i\export PATH="/opt/nvidia/hpc_sdk/Linux_x86_64/2022/compilers/bin/:$PATH"' ~/.bashrc

source ~/.bashrc

# installing healpy and f90wrap

pip3 install healpy
pip3 install f90wrap

```

## Compilation

```
./compile.sh # for the CPU version in docker
```

```
./compile.sh -DGPU=on # for the GPU version in docker
```

```
./compile_ifort.sh # for general CPU version without docker (might need to be modified to incorporate different machine)
```

A known Issue for fastSHT-CPU without docker:

If intel oneapi is installed with a user account, then one may need to run the following command before compiling:
```
export MKL_DIR=~/lib/cmake/mkl-xxxx.x.x/
```
where xxxx.x.x is the mkl version number.

# Examples and Testing

First go to folder ``scripts'', and then:

A comprehensive test and accuracy comparisons with Healpy (may take a long time to run):
```
python test_comprehensive.py
```

A benchmark code:
```
python benchmarks.py
```
or specify the parameters in order of ``nside nsim n_proc niter comparison_flag'':
```
python benchmarks.py 128 1000 8 3 false
```

A test-and-benchmark code for the fix-EB job:
```
python test_fixEB.py
```
or specify the parameters in order of ``nside nsim n_proc niter comparison_flag'':
```
python test_fixEB.py 128 200 8 3 true
```

Notebook that demonstrates the basic interfaces is in  `scripts/demo.ipynb`.

# FAQs

## Linking errors associated with Intel MKL (for installation without docker)

Try pre-load some MKL libraries by

`export LD_PRELOAD=:/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_core.so:/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_lp64.so:/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_thread.so:/opt/intel/oneapi/compiler/latest/linux/compiler/lib/intel64_lin/libiomp5.so`

where ``/opt/intel'' is for the case of installing oneapi with root. If oneapi is installed in a user account, then ``/opt/intel'' can be ``/home/user_name/intel/'' or ``~/intel''

# Citing fastSHT

``Accelerating spherical harmonic transforms for a large number of sky maps'', Chi Tian, Siyu Li, and Hao Liu, https://arxiv.org/abs/2208.10154
