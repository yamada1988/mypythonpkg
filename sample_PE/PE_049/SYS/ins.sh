#!/bin/bash

module load intel
index=49
na=19
gro=PE100.gro

i=`expr $index + $1`
i0=`printf %04d $i`
out=system$i0.gro
L=7.9064
gmx insert-molecules -f $gro -ci $gro -nmol $na -o $out -box $L -radius 0.060 -seed $i -try 20000
