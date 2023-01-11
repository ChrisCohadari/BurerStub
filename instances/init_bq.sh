#!/bin/bash

file="list_filenames_bq.txt"
while read -r line; do 
  [[ "$line" =~ ^#.*$ ]] && continue
  wget --directory-prefix="./bq" "${line}" 
done < "$file"

