#!/bin/bash

cd $PBS_O_WORKDIR
mkdir MD/
n1=0001
gpuindex=0
python ./script/min.py stage0 $n1 $gpuindex  >> OUT/log$n1.txt 

cp -f MD/min$n1.gro SYS/system$n1.gro

