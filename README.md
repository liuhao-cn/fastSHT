fastSHT is a very fast toolkit for doing spherical harmonic transforms on a large number of spherical maps. It converts massive SHT operations to a BLAS level 3 problem and uses the highly optimized matrix multiplication toolkit to accelerate the computation. GPU acceleration is supported and can be very effective. Core code is written in fortran but a Python wrapper is provided and recommended.


# Dependencies

Fortran compiler: `ifort` is recommanded for the CPU version; `nvfortran` is required for the GPU version

Intel MKL library

[`f90wrap`](https://github.com/jameskermode/f90wrap)

`Python3`, `numpy`, `CMake`

# Installation

## (Recommended) Download and compile with the docker image compatable with both CPU and GPU version

```
sudo docker pull rectaflex/intel_nvidia_sdk
```

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

## Compilation

```
./compile.sh # for the CPU version in docker
```

```
./compile.sh -DGPU=on # for the GPU version in docker
```

```
./compile_ifort.sh # for general CPU version (might need to be modified to incorporate different machine)
```

A known Issue:

If intel oneapi is installed with a user account, then one may need to run the following command before compiling:
```
export MKL_DIR=~/lib/cmake/mkl-xxxx.x.x/
```
where xxxx.x.x is the mkl version number.

# Examples and Testing
General tests and comparisons with Healpy is in `scripts/test_all.py`.

Notebook that demonstrates the basic interfaces is in  `scripts/demo.ipynb`.

# FAQs

## Linking errors associated with Intel MKL

Try pre-load some MKL libraries by

`export LD_PRELOAD=:/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_core.so:/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_lp64.so:/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_thread.so:/opt/intel/oneapi/compiler/latest/linux/compiler/lib/intel64_lin/libiomp5.so`

where ``/opt/intel'' is for the case of installing oneapi with root. If oneapi is installed in a user account, then ``/opt/intel'' can be ``/home/user_name/intel/'' or ``~/intel''

# Citing fastSHT

...
