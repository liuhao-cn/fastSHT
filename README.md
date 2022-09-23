# 0. Introduction

fastSHT is a very fast toolkit for doing spherical harmonic transforms on a large number of spherical maps. It converts massive SHT operations to a BLAS level 3 problem and uses the highly optimized matrix multiplication toolkit to accelerate the computation. GPU acceleration is supported and can be very effective. Core code is written in fortran but a Python wrapper is provided and recommended.

To ensure a precise result, fastSHT uses double precision floating numbers (FP64) by default, which prefers a GPU hardware with high FP64 performance. Therefore, the currently best choice is NVIDIA A100 (til Aug-2022), which provides a full-speed FP64 computation with its tensor cores (same performance for double and single precisions).

More details can be found in the following work: 

``Accelerating spherical harmonic transforms for a large number of sky maps'', Chi Tian, Siyu Li, and Hao Liu, https://arxiv.org/abs/2208.10154

# 1. Dependencies

Fortran compiler: `ifort` is recommanded for the CPU version; `nvfortran` is required for the GPU version. The softwares and their versions we used for our tests are listed below:

Intel One API (2022.0.2)

Nvidia HPC SDK (22.3)

`f90wrap (v0.2.7)`  (https://github.com/jameskermode/f90wrap)

`Python3 (3.9.7)`, `numpy (1.21.2)`, `CMake(3.22.1)`

# 2. Installation


## 2.1. Before compilation 

Dependencies and system PATH should be configured properly. 

One can choose either 2.1.1 (with docker) or 2.1.2 (no docker), but there is no need to do both.

### 2.1.1. Environment preparation with docker (much easier)

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
Finally, run the docker image with the following command, which also makes /home available in docker, so if one clone the fastSHT repository to the home directory, it will be available in docker.
```
sudo docker run -it -v /home:/home --gpus all rectaflex/intel_nvidia_sdk
```

### 2.1.2. Environment preparation for ubuntu (no docker)

For non-docker users, we give a sample script for building a compilation environment for fastSHT based on an ubuntu-20.04.

#### 2.1.2.1. Install recent cmake

```
sudo apt update

sudo apt -y install pkg-config sudo build-essential

sudo apt-get install apt-transport-https ca-certificates gnupg software-properties-common wget

wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | sudo apt-key add -

sudo apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main'

sudo apt-get update

sudo apt-get install cmake
```

#### 2.1.2.2. Install Intel oneapi

Choose one way from 2.1.2.2.a and 2.1.2.2.b to Install Intel oneapi. The latter is a partial installation that uses less disk space.

#### 2.1.2.2.a. Install oneapi, way 1: via sudo apt install

```
wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \ | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null

echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list

sudo apt-get update

sudo apt install intel-basekit intel-hpckit

sed -i '1 i\source /opt/intel/oneapi/setvars.sh > /dev/null' ~/.bashrc
source ~/.bashrc
```

#### 2.1.2.2.b. Install oneapi, way 2: a partial installation that saves disk space

```
# Download the Intel oneapi installation packages
wget https://registrationcenter-download.intel.com/akdlm/irc_nas/18673/l_BaseKit_p_2022.2.0.262_offline.sh
wget https://registrationcenter-download.intel.com/akdlm/irc_nas/18679/l_HPCKit_p_2022.2.0.191_offline.sh

# If python 3 is already installed, use the following command to install the basekit (default):
sudo sh l_BaseKit_p_2022.2.0.262_offline.sh -a --silent --eula accept --components intel.oneapi.lin.dpcpp_dbg:intel.oneapi.lin.dpl:intel.oneapi.lin.dpcpp-cpp-compiler:intel.oneapi.lin.mkl.devel:intel.oneapi.lin.tbb.devel

# If python is not installed, use this command to include Intel python in installation
# sudo sh l_BaseKit_p_2022.2.0.262_offline.sh -a --silent --eula accept --components intel.oneapi.lin.dpcpp_dbg:intel.oneapi.lin.dpl:intel.oneapi.lin.dpcpp-cpp-compiler:intel.oneapi.lin.mkl.devel:intel.oneapi.lin.tbb.devel:intel.oneapi.lin.python3

# use the following command to install the HPC kit
sudo sh l_HPCKit_p_2022.2.0.191_offline.sh -a --silent --eula accept --components intel.oneapi.lin.mpi.devel:intel.oneapi.lin.ifort-compiler:intel.oneapi.lin.dpcpp-cpp-compiler-pro

sed -i '1 i\source /opt/intel/oneapi/setvars.sh > /dev/null' ~/.bashrc

source ~/.bashrc
```

#### 2.1.2.3. Install healpy and f90wrap

```
pip3 install healpy f90wrap
```

If no GPU is going to be employed, one can stop here and compile the fastSHT with CPU only using the script `./compile_ifort.sh`. 
If some MKL linking errors are encounted, see FAQs for solutions.

#### 2.1.2.4. Continue with GPU support and install the nvidia hpc sdk:

```
echo 'deb [trusted=yes] https://developer.download.nvidia.com/hpc-sdk/ubuntu/amd64 /' | sudo tee /etc/apt/sources.list.d/nvhpc.list

sudo apt-get update -y

sudo apt-get install -y nvhpc-22-3

sed -i '1 i\export PATH="/opt/nvidia/hpc_sdk/Linux_x86_64/22.3/compilers/bin/:$PATH"' ~/.bashrc

source ~/.bashrc

```

## 2.2. Compilation

If BOTH Intel oneAPI and Nvidia SDK are installed, use

```
./compile.sh # for the CPU version
```
or
```
./compile.sh -DGPU=on # for the GPU version
```

If only Intel oneAPI is installed, use
```
./compile_ifort.sh # for general CPU version (see FAQs if linking errors are encountered)
```

## 2.3 Check the preloads

On some systems one needs to define the preload path; however, on some other systems this should be avoided. Therefore, one should try adding the following preload paths only if one encounters linking errors while importing the SHT python module (mostly these errors regard the openmp libraries). 
```
export LD_PRELOAD=:$LD_PRELOAD:/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_core.so:/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_lp64.so:/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_thread.so:/opt/intel/oneapi/compiler/latest/linux/compiler/lib/intel64_lin/libiomp5.so
```
Here `/opt/intel` is for the case of installing oneapi with root. If oneapi is installed with a user account, then `/opt/intel` should be `/home/user_name/intel/` or `~/intel`.

If, by the above test, this script is found to be necessary, then one should consider adding it to ~/.bashrc for future convenience. You are very much welcome to open an issue if you have experienced other complilation diffuculties. 


# 3. Examples and Testing

First go to folder ``scripts'', and then:

A python file that demonstrates the basic interfaces:
```
python demo.py
```

A comprehensive test and accuracy comparisons with Healpy (may take a long time to run):
```
python test_comprehensive.py
```

A benchmark code:
```
python benchmarks.py  # with default parameters
python benchmarks.py 128 1000 8 3 t2alm false # with parameters in order of nside nsim n_proc niter type comparison_flag
```

A test-and-benchmark code for the fix-EB job:
```
python test_fixEB.py  # with default parameters
python test_fixEB.py 128 200 8 3 true # with parameters in order of nside nsim n_proc niter comparison_flag 
```

# 4. FAQs

## 4.1. Linking errors associated with Intel MKL when import the SHT module (for installation without docker)

### 4.1.a Without GPU
Try pre-load some MKL libraries by

`export LD_PRELOAD=:/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_core.so:/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_lp64.so:/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_thread.so:/opt/intel/oneapi/compiler/latest/linux/compiler/lib/intel64_lin/libiomp5.so`

where `/opt/intel` is for the case of installing oneapi with root. If oneapi is installed in a user account, then `/opt/intel` can be `/home/user_name/intel/` or `~/intel`

### 4.1.b With GPU
Similar to the above, but the last item 

`/opt/intel/oneapi/compiler/latest/linux/compiler/lib/intel64_lin/libiomp5.so` 

should be replaced by 

`/opt/nvidia/hpc_sdk/Linux_x86_64/22.3/REDIST/compilers/lib/libomp.so`

or 

`/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/REDIST/compilers/lib/libomp.so`

where 22.3 or 22.7 depends on the nvidia hpc sdk version.

## 4.2. A known Issue for fastSHT (CPU only) without docker:

If intel oneapi is installed with a user account, then one may need to run the following command before compiling:
```
export MKL_DIR=~/lib/cmake/mkl-xxxx.x.x/
```
where xxxx.x.x is the mkl version number.

# 5. Citing fastSHT

``Accelerating spherical harmonic transforms for a large number of sky maps'', Chi Tian, Siyu Li, and Hao Liu, https://arxiv.org/abs/2208.10154
