#!/bin/bash
#$ -cwd
#$ -l q_node=1
#$ -l h_rt=4:30:00 
#$ -N mdrun
#$ -t 1-1
#$ -o OUT/mdrun.e
#$ -e OUT/mdrun.o


# Module コマンドの初期化
. /etc/profile.d/modules.sh
# CUDA 環境の読込
module load cuda/8.0.61

source ~/.bash_profile

i=$SGE_TASK_ID
num1=$i
num1=`printf %03d $num1`

echo $num1 #$num2 $num3 $num4

export OMP_NUM_THREADS=2

for n1 in $num1 #$num2 $num3 $num4
do

date > OUT/eqlog$n1.txt
starttime=`date "+%s"`
rm OUT/log$n1.txt

done

n1=$num1

###########
# GMX_MIN #
###########
gmx_mpi grompp \
 -f ../../../mdp/min.mdp \
 -c SYS/system$n1.gro \
 -p SYS/topol.top \
 -o MD/gmx$n1.tpr 

rm -f \#*
mpiexec -n 2 gmx_mpi mdrun -deffnm MD/gmx$n1

rm -f MD/\#*
cp MD/gmx$n1.gro SYS/system$n1.gro
python ./script/mdrun.py stage0 $n1 >> OUT/log$n1.txt

    gmx_mpi trjconv \
    -f MD/npt$n1.xtc \
    -o SYS/system$n1.gro \
    -s MD/em$n1.pdb \
    -dump 900 << EOF
0
EOF

rm -f SYS/\#*
cat MD/npt$n1.log >> OUT/density$n1.txt

               
echo "=======" >> OUT/log$n1.txt
echo " stage1 "   >> OUT/log$n1.txt
echo "=======" >> OUT/log$n1.txt

python ./script/mdrun.py stage1 $n1 >> OUT/log$n1.txt

    gmx_mpi trjconv \
    -f MD/npt$n1.xtc \
    -o SYS/system$n1.gro \
    -s MD/em$n1.pdb \
    -dump  1 << EOF
0
EOF
  rm -f SYS/\#*

cat MD/npt$n1.log >> OUT/density$n1.txt

echo "=======" >> OUT/log$n1.txt
echo " stage2 "   >> OUT/log$n1.txt
echo "=======" >> OUT/log$n1.txt

python ./script/mdrun.py stage2 $n1 >> OUT/log$n1.txt

    gmx_mpi trjconv \
    -f MD/npt$n1.xtc \
    -o SYS/system$n1.gro \
    -s MD/em$n1.pdb \
    -dump  1 << EOF
0
EOF
  rm -f SYS/\#*

cat MD/npt$n1.log >> OUT/density$n1.txt

echo "=======" >> OUT/log$n1.txt
echo " stage3 "   >> OUT/log$n1.txt
echo "=======" >> OUT/log$n1.txt

python ./script/mdrun.py stage3 $n1 >> OUT/log$n1.txt

    gmx_mpi trjconv \
    -f MD/npt$n1.xtc \
    -o SYS/system$n1.gro \
    -s MD/em$n1.pdb \
    -dump  50 << EOF
0
EOF
  rm -f SYS/\#*

cat MD/npt$n1.log >> OUT/density$n1.txt

echo "=======" >> OUT/log$n1.txt
echo " stage1 "   >> OUT/log$n1.txt
echo "=======" >> OUT/log$n1.txt

python ./script/mdrun.py stage1 $n1 >> OUT/log$n1.txt

    gmx_mpi trjconv \
    -f MD/npt$n1.xtc \
    -o SYS/system$n1.gro \
    -s MD/em$n1.pdb \
    -dump  1 << EOF
0
EOF
  rm -f SYS/\#*

cat MD/npt$n1.log >> OUT/density$n1.txt

echo "=======" >> OUT/log$n1.txt
echo " stage5 "   >> OUT/log$n1.txt
echo "=======" >> OUT/log$n1.txt

python ./script/mdrun.py stage5 $n1 >> OUT/log$n1.txt

    gmx_mpi trjconv \
    -f MD/npt$n1.xtc \
    -o SYS/system$n1.gro \
    -s MD/em$n1.pdb \
    -dump  1 << EOF
0
EOF
  rm -f SYS/\#*

cat MD/npt$n1.log >> OUT/density$n1.txt
echo "=======" >> OUT/log$n1.txt
echo " stage6 "   >> OUT/log$n1.txt
echo "=======" >> OUT/log$n1.txt

python ./script/mdrun.py stage6 $n1 >> OUT/log$n1.txt

    gmx_mpi trjconv \
    -f MD/npt$n1.xtc \
    -o SYS/system$n1.gro \
    -s MD/em$n1.pdb \
    -dump  50 << EOF
0
EOF
  rm -f SYS/\#*

cat MD/npt$n1.log >> OUT/density$n1.txt

echo "=======" >> OUT/log$n1.txt
echo " stage1 "   >> OUT/log$n1.txt
echo "=======" >> OUT/log$n1.txt

python ./script/mdrun.py stage1 $n1 >> OUT/log$n1.txt

    gmx_mpi trjconv \
    -f MD/npt$n1.xtc \
    -o SYS/system$n1.gro \
    -s MD/em$n1.pdb \
    -dump  1 << EOF
0
EOF
  rm -f SYS/\#*

cat MD/npt$n1.log >> OUT/density$n1.txt

echo "=======" >> OUT/log$n1.txt
echo " stage5 "   >> OUT/log$n1.txt
echo "=======" >> OUT/log$n1.txt

python ./script/mdrun.py stage5 $n1 >> OUT/log$n1.txt

    gmx_mpi trjconv \
    -f MD/npt$n1.xtc \
    -o SYS/system$n1.gro \
    -s MD/em$n1.pdb \
    -dump  1 << EOF
0
EOF
  rm -f SYS/\#*

cat MD/npt$n1.log >> OUT/density$n1.txt
echo "=======" >> OUT/log$n1.txt
echo " stage7 "   >> OUT/log$n1.txt
echo "=======" >> OUT/log$n1.txt

python ./script/mdrun.py stage7 $n1 >> OUT/log$n1.txt

    gmx_mpi trjconv \
    -f MD/npt$n1.xtc \
    -o SYS/system$n1.gro \
    -s MD/em$n1.pdb \
    -dump  1000 << EOF
0
EOF
  rm -f SYS/\#*

cat MD/npt$n1.log >> OUT/density$n1.txt
