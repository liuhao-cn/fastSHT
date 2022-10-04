
# 0. Introduction

fastSHT is a fast toolkit for doing spherical harmonic transforms on a large number of spherical maps. It converts massive SHT operations to a BLAS level 3 problem and uses the highly optimized matrix multiplication toolkit to accelerate the computation. GPU acceleration is supported and can be very effective. The core code is written in Fortran, but a Python wrapper is provided and recommended.

To ensure a precise result, fastSHT uses double precision floating numbers (FP64) by default, which prefers GPU hardware with high FP64 performance. Therefore, the current best choice is the NVIDIA A100 (till Aug-2022), which provides a full-speed FP64 computation with its tensor cores (same performance for double and single precisions).

More technical details can be found in the following work: 

''Accelerating spherical harmonic transforms for a large number of sky maps'', Chi Tian, Siyu Li, and Hao Liu, https://arxiv.org/abs/2208.10154

# 1)  Dependencies

Fortran compiler `ifort` is recommanded for the CPU version, and `nvfortran` is required for the GPU version. The main dependencies and their versions are listed below:

`Intel One API (2022.0.2)`

`Nvidia HPC SDK (22.3 or 22.7)`

`f90wrap (v0.2.7)`  (https://github.com/jameskermode/f90wrap)

`Python3 (3.9.7)`

`numpy (1.21.2)`

`CMake(3.22.1)`

# 2) Environment configuration

The build dependencies and system PATH should be configured properly.  This can be done by following either 2.1 (with docker) or 2.2 (no docker).

First of all, clone the repository:
```
git clone https://github.com/liuhao-cn/fastSHT.git
```

## 2.1) Environment configuration with docker

Here we assume docker is already installed and available. See https://docs.docker.com/engine/install/ for a docker installation instruction.

The first step of configuration with docker is to prepare the docker image. This can be done in one of the following two ways:

A) Build the image locally from a Dockerfile (uses less disk space):
```
cd ./docker
sudo docker build -t fastsht:gpu . 
```
B) Pull the pre-built image (uses more disk space):
```
sudo docker pull rectaflex/intel_nvidia_sdk
sudo docker image tag rectaflex/intel_nvidia_sdk fastsht:gpu
```

When the docker image is prepared, one needs to install the Nvidia container runtime by

```
curl -s -L https://nvidia.github.io/nvidia-container-runtime/gpgkey |   sudo apt-key add -

distribution=$(. /etc/os-release;echo $ID$VERSION_ID)

curl -s -L https://nvidia.github.io/nvidia-container-runtime/$distribution/nvidia-container-runtime.list |   sudo tee /etc/apt/sources.list.d/nvidia-container-runtime.list

sudo apt-get update

sudo apt-get install nvidia-container-runtime

sudo systemctl restart docker
```
Finally, run the docker image with the following command: 
```
sudo docker run -it -v /home:/home --gpus all fastsht:gpu
```
which makes all GPUs available in docker and also makes `/home` available, so if one clones the fastSHT repository to `/home` on the host machine, it will be available in docker.

## 2.2) Environment configuration without docker (for ubuntu 20.04)

The environment configuration without docker is tested for ubuntu 20.04, and can be done either with the auto-configuration script or a step-by-step manual configuration.

### 2.2.1) Auto-configuration script

To automatically install necessary environment, run 
```
./configure.sh
``` 
The default behavior of `configure.sh` is a FULL installation of Intel ONE API and NO CUDA.  Use the argument '--part-no-python' for part installtaion of ONE API with no intel python, and '--part-with-python' for part install with the intel python. Use the argument '--with-gpu' to install NVIDIA HPC SDK to enable GPU capability.

### 2.2.2) Step-by-step manual configuration

#### step-1) Install recent cmake
```
sudo apt update

sudo apt -y install pkg-config sudo build-essential

sudo apt-get install apt-transport-https ca-certificates gnupg software-properties-common wget

wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | sudo apt-key add -

sudo apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main'

sudo apt-get update

sudo apt-get install cmake
```

#### step-2) Install Intel oneapi

Choose one of the following two ways to Install Intel oneapi. Way 2 is a partial installation that uses less disk space.

##### step-2-way 1: via sudo apt install

```
wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \ | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null

echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list

sudo apt-get update

sudo apt install intel-basekit intel-hpckit

sed -i '1 i\source /opt/intel/oneapi/setvars.sh > /dev/null' ~/.bashrc

source ~/.bashrc
```

##### step-2-way 2: a partial installation that saves disk space

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

#### step-3) Install healpy and f90wrap

```
pip3 install healpy f90wrap
```

If no GPU is going to be employed, one can stop here and jump to compilation.
 
If some MKL linking errors are encounted, see FAQs for solutions.

#### step-4) Continue with GPU support and install the nvidia hpc sdk:

```
echo 'deb [trusted=yes] https://developer.download.nvidia.com/hpc-sdk/ubuntu/amd64 /' | sudo tee /etc/apt/sources.list.d/nvhpc.list

sudo apt-get update -y

sudo apt-get install -y nvhpc-22-3

sed -i '1 i\export PATH="/opt/nvidia/hpc_sdk/Linux_x86_64/22.3/compilers/bin/:$PATH"' ~/.bashrc

source ~/.bashrc
```

# 3) Compilation

Use the following command to compile the CPU version:
```
./compile.sh
```
or the following command to compile the GPU version:
```
./compile.sh -DGPU=on
```


# 4) Examples and test scripts

First `cd ./scripts` and then run one of the following scripts: 

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

# 5) FAQs

## 5.1) Linking errors associated with Intel MKL

For a non-docker installation, a linking error may occur when attempting to import the python module `SHT`. This can be fixed by defining the preload paths:

### The CPU-only case
Try pre-load some MKL libraries by

```
export LD_PRELOAD=:/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_core.so:/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_lp64.so:/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_thread.so:/opt/intel/oneapi/compiler/latest/linux/compiler/lib/intel64_lin/libiomp5.so
```

where `/opt/intel` is for the case of installing oneapi with root. If oneapi is installed in a user account, then `/opt/intel` can be `/home/user_name/intel/` or `~/intel`

### The case with GPU

Similar to the above, but the last item 

```
/opt/intel/oneapi/compiler/latest/linux/compiler/lib/intel64_lin/libiomp5.so
``` 

should be replaced by 

```
/opt/nvidia/hpc_sdk/Linux_x86_64/22.3/REDIST/compilers/lib/libomp.so
```

or 

`/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/REDIST/compilers/lib/libomp.so`

where 22.3 or 22.7 depends on the nvidia hpc sdk version. Because the GPU version uses a different omp library.

If, the above preload script is found to be necessary, then one should consider adding it to ~/.bashrc for future convenience. You are very much welcome to open an issue if you have experienced other complilation diffuculties. 

## 5.2. A known Issue for fastSHT (CPU only) without docker:

If intel oneapi is installed with a user account, then one may need to run the following command before compiling:
```
export MKL_DIR=~/lib/cmake/mkl-xxxx.x.x/
```
where xxxx.x.x is the mkl version number.

# 6. Citing fastSHT

``Accelerating spherical harmonic transforms for a large number of sky maps'', Chi Tian, Siyu Li, and Hao Liu, https://arxiv.org/abs/2208.10154
