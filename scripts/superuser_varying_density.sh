#!/bin/bash

# default parameters
querysize=9
windowsize=30000
density=0.50

resultdir=../results/superuser_varying_density
tcm=../TCM_snap

rm -rf ${resultdir}
mkdir -p ${resultdir}

data=../datasets/superuser.graph
querydir=../querysets/superuser

for d in 0.00 0.25 0.50 0.75 1.00
do
  for n in {0..99}
  do 
    q=${querysize}_${n}_${d}
    timeout 1h ${tcm} ${data} ${querydir}/${q} ${windowsize} > ${resultdir}/${q}_${windowsize}
  done
done
