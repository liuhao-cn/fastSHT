fastSHT is a very fast toolkit for doing spherical harmonic transforms on a large number of spherical maps. It converts massive SHT operations to a BLAS level 3 problem and uses the highly optimized matrix multiplication toolkit to accelerate the computation. GPU acceleration is supported and can be very effective. Core code is written in fortran but a Python wrapper is provided and recommended.


# Dependencies

Fortran compiler: `ifort` is recommanded for the CPU version; `nvfortran` is required for the GPU version

Intel MKL library

[`f90wrap`](https://github.com/jameskermode/f90wrap)

`Python3`, `numpy`, `CMake`


# (Recommended) Download and compile with the docker image that is compatable with both CPU and GPU version

```
docker pull rectaflex/intel_nvidia_sdk
```

# Large file support
Some data file needs to be downloaded using git-lfs. Please refer to the git-lfs manual for more details. 

For example: in Ubuntu this can be done by
```
sudo apt install git-lfs
git pull
git-lfs pull
```

# Compilation

```
./compile.sh # for the CPU version
```

```
./compile.sh -DGPU=on # for the GPU version
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

`export LD_PRELOAD=:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_core.so:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_intel_lp64.so:/opt/intel/oneapi/mkl/2022.0.2/lib/intel64/libmkl_intel_thread.so:/opt/intel/oneapi/compiler/2022.0.2/linux/compiler/lib/intel64_lin/libiomp5.so`

where ``/opt/intel'' is for the case of installing oneapi with root. If oneapi is installed in a user account, then this can be /home/user_name/intel/

# Citing fastSHT

...
