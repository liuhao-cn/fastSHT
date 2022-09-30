#!/bin/bash

p1=all
p2=no

case $1 in
    --part-no-python)
    p1=--part-no-python
    ;;
    --part-with-python)
    p1=--part-with-python
    ;;
    --with-gpu)
    p2=--with-gpu
    ;;
esac

case $2 in
    --part-no-python)
    p1=--part-no-python
    ;;
    --part-with-python)
    p1=--part-with-python
    ;;
    --with-gpu)
    p2=--with-gpu
    ;;
esac


PKG_OK=$(which ifort|grep "ifort")

if [[ "$PKG_OK" == "" ]]; then

case $p1 in
    --part-no-python)
	echo "Full install intel-one-api with no python"
	# This script will only install the necessary parts of Intel oneapi, including Intel python.

	# Install the most recent cmake and build tools
	sudo apt update

	sudo apt -y install gfortran gcc python3-pip pkg-config sudo build-essential

	sudo apt-get -y install apt-transport-https ca-certificates gnupg software-properties-common wget


	wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | sudo apt-key add -
	sudo apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main'
	sudo apt-get update
	sudo apt-get -y install cmake


	# Download the Intel oneapi installation packages
	wget https://registrationcenter-download.intel.com/akdlm/irc_nas/18673/l_BaseKit_p_2022.2.0.262_offline.sh -c

	wget https://registrationcenter-download.intel.com/akdlm/irc_nas/18679/l_HPCKit_p_2022.2.0.191_offline.sh -c


	# If python 3 is already installed, use the following command to install the basekit without python:
	sudo sh l_BaseKit_p_2022.2.0.262_offline.sh -a --silent --eula accept --components intel.oneapi.lin.dpcpp_dbg:intel.oneapi.lin.dpl:intel.oneapi.lin.dpcpp-cpp-compiler:intel.oneapi.lin.mkl.devel:intel.oneapi.lin.tbb.devel

	# # If python is not installed, use this command to include Intel python in installation
	# sudo sh l_BaseKit_p_2022.2.0.262_offline.sh -a --silent --eula accept --components intel.oneapi.lin.dpcpp_dbg:intel.oneapi.lin.dpl:intel.oneapi.lin.dpcpp-cpp-compiler:intel.oneapi.lin.mkl.devel:intel.oneapi.lin.tbb.devel:intel.oneapi.lin.python3

	# use the following command to install the HPC kit
	sudo sh l_HPCKit_p_2022.2.0.191_offline.sh -a --silent --eula accept --components intel.oneapi.lin.mpi.devel:intel.oneapi.lin.ifort-compiler:intel.oneapi.lin.dpcpp-cpp-compiler-pro

	sed -i '1 i\source /opt/intel/oneapi/setvars.sh > /dev/null' ~/.bashrc

	source ~/.bashrc


	# install healpy and f90wrap
	pip3 install healpy f90wrap numba

	source ~/.bashrc

        ;;
    --part-with-python)
	echo "Full install intel-one-api with intel python"
	# This script will only install the necessary parts of Intel oneapi, including Intel python.

	# Install the most recent cmake and build tools
	sudo apt update

	sudo apt -y install gfortran gcc python3-pip pkg-config sudo build-essential

	sudo apt-get -y install apt-transport-https ca-certificates gnupg software-properties-common wget

	wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | sudo apt-key add -
	sudo apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main'
	sudo apt-get update
	sudo apt-get -y install cmake


	# Download the Intel oneapi installation packages
	wget https://registrationcenter-download.intel.com/akdlm/irc_nas/18673/l_BaseKit_p_2022.2.0.262_offline.sh -c

	wget https://registrationcenter-download.intel.com/akdlm/irc_nas/18679/l_HPCKit_p_2022.2.0.191_offline.sh -c


	# If python 3 is already installed, use the following command to install the basekit without python:
	# sudo sh l_BaseKit_p_2022.2.0.262_offline.sh -a --silent --eula accept --components intel.oneapi.lin.dpcpp_dbg:intel.oneapi.lin.dpl:intel.oneapi.lin.dpcpp-cpp-compiler:intel.oneapi.lin.mkl.devel:intel.oneapi.lin.tbb.devel

	# # If python is not installed, use this command to include Intel python in installation
	sudo sh l_BaseKit_p_2022.2.0.262_offline.sh -a --silent --eula accept --components intel.oneapi.lin.dpcpp_dbg:intel.oneapi.lin.dpl:intel.oneapi.lin.dpcpp-cpp-compiler:intel.oneapi.lin.mkl.devel:intel.oneapi.lin.tbb.devel:intel.oneapi.lin.python3

	# use the following command to install the HPC kit
	sudo sh l_HPCKit_p_2022.2.0.191_offline.sh -a --silent --eula accept --components intel.oneapi.lin.mpi.devel:intel.oneapi.lin.ifort-compiler:intel.oneapi.lin.dpcpp-cpp-compiler-pro

	sed -i '1 i\source /opt/intel/oneapi/setvars.sh > /dev/null' ~/.bashrc

	source ~/.bashrc


	# install healpy and f90wrap
	pip3 install healpy f90wrap numba

	source ~/.bashrc

        shift
        ;;
    *)
	echo "Full install intel-one-api"
	sudo apt update

	sudo apt -y install gfortran

	sudo apt -y install pkg-config sudo build-essential

	sudo apt-get -y install apt-transport-https ca-certificates gnupg software-properties-common wget


	wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | sudo apt-key add -
	sudo apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main'
	sudo apt-get update
	sudo apt-get -y install cmake


	wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \ | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null

	echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list

	sudo apt-get update

	sudo apt -y install intel-basekit intel-hpckit

	sed -i '1 i\source /opt/intel/oneapi/setvars.sh > /dev/null' ~/.bashrc

	source ~/.bashrc

	pip3 install healpy f90wrap

	source ~/.bashrc
        ;;


esac

fi

PKG_OK=$(which ifort|grep "ifort")

if [[ "$PKG_OK" == "" ]]; then

    echo "=====================Intel One API installation failed================"
    exit
fi



case $p2 in
    --with-gpu)
	echo 'deb [trusted=yes] https://developer.download.nvidia.com/hpc-sdk/ubuntu/amd64 /' | sudo tee /etc/apt/sources.list.d/nvhpc.list

	sudo apt-get update -y

	sudo apt-get install -y nvhpc-22-3

	sed -i '1 i\export PATH="/opt/nvidia/hpc_sdk/Linux_x86_64/22.3/compilers/bin/:$PATH"' ~/.bashrc

	source ~/.bashrc
	PKG_OK=$(which nvfortran|grep "nvfortran")

	if [[ "$PKG_OK" == "" ]]; then

	    echo "=====================NV HPC SDK installation failed================"
	    exit
	fi

	
	if [[ "$p1" == "--part-with-python" ||  "$p1" == "--part-no-python" ]]; then
	       sed -i '1 i\export LD_PRELOAD=:/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_core.so:/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_lp64.so:/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_thread.so:/opt/nvidia/hpc_sdk/Linux_x86_64/22.3/REDIST/compilers/lib/libomp.so' ~/.bashrc
	fi
    
	source ~/.bashrc
	;;
	 
    *)
	echo "no NV HPC SDK is going to install"

	if [[ "$p1" == "--part-with-python" ||  "$p1" == "--part-no-python" ]]; then
	       sed -i '1 i\export LD_PRELOAD=:/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_core.so:/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_lp64.so:/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_thread.so:/opt/intel/oneapi/compiler/latest/linux/compiler/lib/intel64_lin/libiomp5.so' ~/.bashrc
	fi

    
esac


    
