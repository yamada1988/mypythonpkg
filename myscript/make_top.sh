inpf='SYS/topol.top'
outf='SYS/processed.top'
inpgro='SYS/system001.gro'
inpmdp='../../../mdp/min.mdp'

gmx_mpi_d grompp -f $inpmdp -c $inpgro -p $inpf -pp $outf

rm -f \#* SYS/\#*

python script/make_top.py $outf

nrep=8
# "effective" temperature range
tmin=453
tmax=753

# build geometric progression
list=$(
awk -v n=$nrep \
    -v tmin=$tmin \
    -v tmax=$tmax \
    'BEGIN{for(i=0;i<n;i++){
    t=tmin*exp(i*log(tmax/tmin)/(n-1));
    printf(t); if(i<n-1)printf(",");
    }
    }'
)
    
echo $list
echo $list > SYS/temperature.dat
for((i=0;i<nrep;i++))
do

  # choose lambda as T[0]/T[i]
  # remember that high temperature is equivalent to low lambda
  lambda=$(echo $list | awk 'BEGIN{FS=",";}{print $1/$'$((i+1))';}')
  # process topology
  j=`expr $i + 1`
  j=`printf %02d $j`
  plumed partial_tempering $lambda < $outf > SYS/topol_$j.top
done
