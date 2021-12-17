#!/bin/bash
#PJM -N "rg16P40replicanum"
#PJM -L "rscunit=ito-b"
#PJM -L "rscgrp=hp-g-1-dbg"
#PJM -L "vnode=1"
#PJM -L "vnode-core=9"
#PJM -L "elapse=00:60:00"
#PJM -X
#PJM -e "OUT/rerun_gromacs.e.replicanum"
#PJM -o "OUT/rerun_gromacs.o.replicanum"

#---- Program Execution --------#                          
#module load cuda/8.0
#module load intel/4017

export OMP_NUM_THREADS=2

i=replicanum
m=4
n1=`expr $m \* $i - $m + 1`
n2=`expr $n1 + $m - 1`

module load gromacs/2020.6-gpu

source ~/.bashrc
mkdir DAT

dir=gromacs/

# grompp
if [ $i == 1 ]; then

  mkdir $dir
  for ind in "B0"
  do
    if [ ! -f $dir/$ind.tpr ] ; then
      mpirun -np 1 gmx_mpi_d grompp -f ../../mdp/md.mdp -c SYS/system0001.gro -p SYS/topol_$ind.top -o $dir/$ind -maxwarn 1
      rm -f *mdout*
      rm -f $dir/\#*
    fi
  done

  for ind in "S0" "S0B0"
  do
    if [ ! -f $dir/$ind.tpr ]; then
      mpirun -np 1 gmx_mpi_d grompp -f ../../mdp/md.mdp -c SYS/solute.gro -p SYS/topol_${ind}slt.top -o $dir/$ind -maxwarn 1
      rm -f *mdout*
      rm -f $dir/\#*
    fi
  done
fi
date

#trjcat
for n in `seq $n1 $n2`
do

  num=`printf %04d $n`

  if [ ! -f MD/slt$num.xtc ]; then
  mpirun -np 1 gmx_mpi_d trjconv -f MD/npt$num.xtc -s gromacs/S0.tpr -o MD/slt$num.xtc << EOF
0
EOF
fi

done

rm -f $dir/\#*

#rerun
for n in `seq $n1 $n2`
do

  for ind in "S0" "S0B0" 
  do
    num=`printf %04d $n`
    echo $n
    mpirun -np 2 gmx_mpi mdrun -s $dir/$ind.tpr -rerun MD/slt$num.xtc -deffnm $dir/npt${num}_$ind 
    rm mdout.mdp
    rm $dir/npt${num}_$ind.trr
  
    mpirun -np 1 gmx_mpi energy -f $dir/npt${num}_$ind.edr -o $dir/npt${num}_$ind.xvg << EOF
9 0
EOF

rm -f $dir/\#*

    sed -e '1,24d' $dir/npt${num}_$ind.xvg > $dir/npt${num}_$ind.dat
    rm $dir/npt${num}_$ind.xvg 
  
  done

  for ind in "B0" 
  do
    num=`printf %04d $n`
    echo $n
    mpirun -np 2 gmx_mpi mdrun -s $dir/$ind.tpr -rerun MD/npt$num.xtc -deffnm $dir/npt${num}_$ind 
    rm mdout.mdp
    rm $dir/npt${num}_$ind.trr

    mpirun -np 1 gmx_mpi energy -f $dir/npt${num}_$ind.edr -o $dir/npt${num}_$ind.xvg << EOF
10 0
EOF

rm -f $dir/\#*

    sed -e '1,24d' $dir/npt${num}_$ind.xvg > $dir/npt${num}_$ind.dat
    rm $dir/npt${num}_$ind.xvg

  done
  


done


rm -f \#* *mdout*

date
