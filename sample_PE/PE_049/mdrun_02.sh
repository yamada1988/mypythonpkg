#!/bin/bash
module load cuda/10.1.2

cd $PBS_O_WORKDIR
n1=0001

gpuindex=0
python ./script/sample.py stage2 $n1 $gpuindex  >> OUT/log$n1.txt 

cp -f MD/npt$n1.gro SYS/system$n1.gro
