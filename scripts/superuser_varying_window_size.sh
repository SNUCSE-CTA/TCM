#!/bin/bash

# default parameters
querysize=9
windowsize=30000
density=0.50

resultdir=../results/superuser_varying_window_size
tcm=../TCM_snap

rm -rf ${resultdir}
mkdir -p ${resultdir}

data=../datasets/superuser.graph
querydir=../querysets/superuser

for wsize in 10000 20000 30000 40000 50000
do
  for n in {0..99}
  do 
    q=${querysize}_${n}_${density}
    timeout 1h ${tcm} ${data} ${querydir}/${q} ${wsize} > ${resultdir}/${q}_${wsize}
  done
done
