#!/bin/bash
#PJM -N "rgB16A034replicanum"
#PJM -L "rscunit=ito-b"
#PJM -L "rscgrp=hp-g-4"
#PJM -L "vnode=1"
#PJM -L "vnode-core=36"
#PJM -L "elapse=02:00:00"
#PJM -X
#PJM -e "OUT/reruns.e.replicanum"
#PJM -o "OUT/reruns.o.replicanum"

#---- Program Execution --------#                          
module load cuda/8.0
module load intel/2017

i=replicanum
m=10
n1=`expr $m \* $i - $m + 1`
n2=`expr $n1 + $m - 1`


source ~/.bashrc
mkdir DAT


date

gpuindex=0
for n in `seq $n1 $n2`
do

echo $n
num=`printf %04d $n`
for ind in "_ghost" "_B0" 
do

  g=`echo $(( $gpuindex % 4 ))`
  ~/software/anaconda3/bin/python script/rerun$ind.py check $num $g &
  echo $ind $gpuindex $g
  gpuindex=`expr $gpuindex + 1`
done

done
wait

date


for n in `seq $n1 $n2`
do

echo $n
~/software/anaconda3/bin/python script/make_SysWght.py $n

done
