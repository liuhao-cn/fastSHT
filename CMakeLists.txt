project(fastSHT)
cmake_minimum_required(VERSION 3.16.2)

option(GPU "GPU" OFF)
option(INTEL "INTEL" ON)

enable_language (Fortran)

if(INTEL AND NOT GPU)
set(CMAKE_Fortran_COMPILER  ifort)
else()
set(CMAKE_Fortran_COMPILER  nvfortran)
endif()

find_package(MKL CONFIG REQUIRED)


file(GLOB fastSHT_SOURCES src/alm2cl_TEB.f90 src/cl2alm_TEB.f90            src/compute_fft.f90           src/convert_alm.f90           src/fft_mapping.f90
 src/filt_alm.f90              src/get_plm_recursive.f90                   src/norm_alm.f90             src/sht_data_module.f90        src/sht_core.f90 src/sht_main.f90 src/sht_data_init.f90)



set(MKLROOT /opt/intel/oneapi/mkl/2022.0.2/)


if(GPU)
message("Will compile with GPU (supported by openacc)")
set(CMAKE_Fortran_FLAGS "-fPIC -O3 -fopenmp -Mcuda -Mcudalib=cublas -lcufft -acc -ta=tesla -fast -cpp -DGPU=TRUE")
else()
message("Will compile with CPU only")
set(CMAKE_Fortran_FLAGS "-fPIC -O3 -fopenmp -cpp")
endif()

add_library(fastSHT_fort STATIC ${fastSHT_SOURCES})

#set(CMAKE_VERBOSE_MAKEFILE ON)

target_link_libraries(fastSHT_fort MKL::MKL )
