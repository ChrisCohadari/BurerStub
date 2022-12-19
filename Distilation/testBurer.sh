#!/bin/bash

instance_directory=~/Documents/studies/Sorbonne/TER/BurerStub/instances
for entry in "${instance_directory}"/*
do 
  ./burer "$entry"
done
