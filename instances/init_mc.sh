#!/bin/bash

file="list_filenames_mc.txt"
while read -r line; do 
  [[ "$line" =~ ^#.*$ ]] && continue
  wget --directory-prefix="./mc" "${line}" 
done < "$file"

