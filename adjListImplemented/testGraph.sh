#!/bin/bash

instance_directory=~/Documents/studies/Sorbonne/TER/BurerStub/instances/mc
for entry in "${instance_directory}"/*
do 
  ./graph "$entry"
done
