#!/bin/bash

instance_directory=../../instances/mc
for Mval in 05 10 15 20 30 40 50
do
  for Nval in 05 10 15
  do 
    echo M=$Mval,N=$Nval >> log.M$Mval.N$Nval.txt
    for entry in "${instance_directory}"/*
    do 
      ./burer.M$Mval.N$Nval "$entry" >> log.M$Mval.N$Nval.txt
    done
  done
done

