#!/bin/bash

cp benchmarks.txt data.txt

for Mval in 05 10 15 20 30 40 50
do
  for Nval in 05 10 15
  do 
      paste data.txt log.M$Mval.N$Nval.txt >> data.txt
      echo Ho
  done
done


