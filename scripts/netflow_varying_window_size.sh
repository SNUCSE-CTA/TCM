#!/bin/bash

# default parameters
querysize=9
windowsize=30000
density=0.50

resultdir=../results/netflow_varying_window_size
tcm=../TCM_network

rm -rf ${resultdir}
mkdir -p ${resultdir}

data=../datasets/netflow.graph
querydir=../querysets/netflow

for wsize in 10000 20000 30000 40000 50000
do
  for n in {0..99}
  do 
    q=${querysize}_${n}_${density}
    timeout 1h ${tcm} ${data} ${querydir}/${q} ${wsize} > ${resultdir}/${q}_${wsize}
  done
done
