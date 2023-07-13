#!/bin/bash

# default parameters
querysize=9
windowsize=30000
density=0.50

D=stackoverflow

resultdir=../results/${D}
tcm=../TCM

rm -rf ${resultdir}
mkdir -p ${resultdir}/querysize
mkdir -p ${resultdir}/density
mkdir -p ${resultdir}/windowsize

data=../datasets/${D}.graph
querydir=../querysets/${D}

for qsize in 5 7 9 11 13 15
do
  for n in {0..99}
  do 
    q=${qsize}_${n}_${density}
    timeout 1h ${tcm} ${data} ${querydir}/${q} ${windowsize} > ${resultdir}/querysize/${q}_${windowsize}
  done
done

for d in 0.00 0.25 0.50 0.75 1.00
do
  for n in {0..99}
  do 
    q=${querysize}_${n}_${d}
    timeout 1h ${tcm} ${data} ${querydir}/${q} ${windowsize} > ${resultdir}/density/${q}_${windowsize}
  done
done

for wsize in 10000 20000 30000 40000 50000
do
  for n in {0..99}
  do 
    q=${querysize}_${n}_${density}
    timeout 1h ${tcm} ${data} ${querydir}/${q} ${wsize} > ${resultdir}/windowsize/${q}_${wsize}
  done
done
