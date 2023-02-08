#!/bin/bash

instance_directory=../instances/mc
for entry in "${instance_directory}"/*
do 
  ./burer.M10.N10 "$entry" >> log.M10.N10.txt
done
