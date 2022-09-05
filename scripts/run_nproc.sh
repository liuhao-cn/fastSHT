#!/bin/bash  

# syntax: 
# ./run_nproc.sh nside nsim nproc_start nproc_end nproc_step
for(( i=$3; i<=$4; i=i+$5 ));  do   
	python nproc.py $1 $2 $i 3
done  
