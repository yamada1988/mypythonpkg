#!/bin/bash

for i in `seq -f %04g 1 10`
do
  cp mixture.inp mixture.inp.bak
  sed -i -e "s/index/$i/g" mixture.inp
  packmol < mixture.inp
  mv mixture.inp.bak mixture.inp
done
