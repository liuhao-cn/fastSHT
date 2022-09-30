FROM ubuntu:20.04

SHELL ["/bin/bash", "-c"]

####################################
# First fix the time zone
####################################
RUN ln -fs /usr/share/zoneinfo/Asia/Hong_Kong /etc/localtime \
    && apt update \
    && apt install -y tzdata \
    && dpkg-reconfigure --frontend noninteractive tzdata \
####################################
# install necessary packages
####################################
    && apt -y install git nano python3 python3-pip apt-utils gfortran gcc pkg-config sudo build-essential apt-transport-https ca-certificates gnupg software-properties-common wget \
####################################
# install recent cmake
####################################
    && wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | apt-key add - \
    && apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main' \
    && apt update \
    && apt -y install cmake \
####################################
# install oneapi (only necessary packages)
####################################
    && wget https://registrationcenter-download.intel.com/akdlm/irc_nas/18673/l_BaseKit_p_2022.2.0.262_offline.sh \
    && chmod 755 l_BaseKit_p_2022.2.0.262_offline.sh \
    && ./l_BaseKit_p_2022.2.0.262_offline.sh -a --silent --eula accept --components intel.oneapi.lin.dpcpp_dbg:intel.oneapi.lin.dpl:intel.oneapi.lin.dpcpp-cpp-compiler:intel.oneapi.lin.mkl.devel:intel.oneapi.lin.tbb.devel\
    && rm l_BaseKit_p_2022.2.0.262_offline.sh \
    && wget https://registrationcenter-download.intel.com/akdlm/irc_nas/18679/l_HPCKit_p_2022.2.0.191_offline.sh \
    && chmod 755 l_HPCKit_p_2022.2.0.191_offline.sh \
    && ./l_HPCKit_p_2022.2.0.191_offline.sh -a --silent --eula accept --components intel.oneapi.lin.mpi.devel:intel.oneapi.lin.ifort-compiler:intel.oneapi.lin.dpcpp-cpp-compiler-pro \
    && rm l_HPCKit_p_2022.2.0.191_offline.sh \
    && sed -i '1 i\source /opt/intel/oneapi/setvars.sh > /dev/null' ~/.bashrc \
####################################
# install nvidia sdk
####################################
    && wget https://developer.download.nvidia.com/hpc-sdk/22.7/nvhpc_2022_227_Linux_x86_64_cuda_11.7.tar.gz \
    && tar xpzf nvhpc_2022_227_Linux_x86_64_cuda_11.7.tar.gz \
    && rm nvhpc_2022_227_Linux_x86_64_cuda_11.7.tar.gz \
    && nvhpc_2022_227_Linux_x86_64_cuda_11.7/install \
    && rm -rf nvhpc_2022_227_Linux_x86_64_cuda_11.7/ \
    && sed -i '1 i\export PATH="/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/compilers/bin/:$PATH"' ~/.bashrc \
    && sed -i '1 i\export LD_PRELOAD=:/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_core.so:/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_lp64.so:/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_thread.so:/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/REDIST/compilers/lib/libomp.so' ~/.bashrc \
    && source ~/.bashrc \
####################################
# install python packages
####################################
    && pip3 install -i https://pypi.tuna.tsinghua.edu.cn/simple healpy f90wrap numba
