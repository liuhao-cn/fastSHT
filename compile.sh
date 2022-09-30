#!/bin/bash


echo ${1}
rm -rf obj/

mkdir obj

cp -r src/ obj/

cd obj


if [[ ${1} == *"GPU"* ]]; then

    export FC=nvfortran

    cmake .. ${1}

    make -j4

    cp *.mod src/

    gfortran -E -cpp -DGPU=True ../src/sht_data_module.f90 > src/sht_data_module.f90
    gfortran -E -cpp -DGPU=True ../src/sht_main.f90 > src/sht_main.f90

    f90wrap -m fastSHT -k src/kind_map  src/sht_main.f90 src/convert_alm.f90 src/sht_data_init.f90

    LDFLAGS="-cpp -Mcuda -acc -Mcudalib=cublas -fPIC -lcufft -pgf90libs -mp -lpthread -lm -ldl -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread" f2py-f90wrap --fcompiler=nv -c -m _fastSHT f90wrap_*.f90 --f90flags="-cpp -acc -Mcud\
a -Mcudalib=cublas -pgf90libs -mp -lpthread -lm -ldl" *.a

else

    cmake .. ${1}

    make -j4

    cp *.mod src/


    gfortran -E -cpp ../src/sht_data_module.f90 > src/sht_data_module.f90
    gfortran -E -cpp  ../src/sht_main.f90 > src/sht_main.f90

    f90wrap -m fastSHT -k src/kind_map  src/sht_main.f90 src/sht_data_init.f90 src/convert_alm.f90

    LDFLAGS="-cpp -lpthread -ldl -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -qopenmp" f2py-f90wrap --fcompiler=intelem -c -m _fastSHT f90wrap_*.f90 --f90flags="-cpp   -lpthread -ldl" *.a
fi
