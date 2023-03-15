#!/bin/bash

cp benchmarks.txt data0.txt

i=0
for Mval in 05 10 15 20 30 40 50
do
  for Nval in 05 10 15
  do 
	paste --delimiters=, data$i.txt logs/log.M$Mval.N$Nval.txt > data$((i+1)).txt
  rm data$i.txt
  ((i++))
  done
done

mv data$i.txt data.txt



