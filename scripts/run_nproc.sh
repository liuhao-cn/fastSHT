#!/bin/bash
niter=3

# test 1
nside=64
nsim=40000
for(( i=8; i<=48; i=i+4 ));  do
	python nproc.py $nside $nsim $i $niter
done

# test 2
nside=128
nsim=20000
for(( i=8; i<=48; i=i+4 ));  do
    python nproc.py $nside $nsim $i $niter
done

# test 3
nside=256
nsim=8000
for(( i=8; i<=48; i=i+4 ));  do
    python nproc.py $nside $nsim $i $niter
done

# test 4
nside=256
nsim=16000
for(( i=8; i<=48; i=i+4 ));  do
    python nproc.py $nside $nsim $i $niter
done

# test 5
nside=512
nsim=4000
for(( i=8; i<=48; i=i+4 ));  do
    python nproc.py $nside $nsim $i $niter
done

# test 6
nside=1024
nsim=2000
for(( i=8; i<=48; i=i+4 ));  do
    python nproc.py $nside $nsim $i $niter
done
