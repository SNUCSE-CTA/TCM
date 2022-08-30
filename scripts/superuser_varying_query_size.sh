#!/bin/bash

# default parameters
querysize=9
windowsize=30000
density=0.50

resultdir=../results/superuser_varying_query_size
tcm=../TCM_snap

rm -rf ${resultdir}
mkdir -p ${resultdir}

data=../datasets/superuser.graph
querydir=../querysets/superuser

for qsize in 5 7 9 11 13 15
do
  for n in {0..99}
  do 
    q=${qsize}_${n}_${density}
    timeout 1h ${tcm} ${data} ${querydir}/${q} ${windowsize} > ${resultdir}/${q}_${windowsize}
  done
done
