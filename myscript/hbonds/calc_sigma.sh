#!/bin/bash


gmx_mpi energy -f MD/npt_rsigma.edr -o sigma.xvg << EOF
32 33 36 0
EOF

rm -f \#*

num=`wc -l sigma.xvg | awk  '{print $1}'`
num=`expr $num - 26`
cp script/calc_sigma.f90 script/run.f90
sed -i -e "s/num/$num/g" script/run.f90

gfortran script/run.f90 -o calc_sigma -fopenmp
./calc_sigma
