#!/bin/bash

cp benchmarks.txt data1.txt

i=0
for Mval in 05 10 15 20 30 40 50
do
  for Nval in 05 10 15
  do 
      if ((i %2 == 0))
      then
	paste --delimiters=, data1.txt log.M$Mval.N$Nval.txt > data2.txt
      fi

      if ((i % 2 == 1))
      then
	paste --delimiters=, data2.txt log.M$Mval.N$Nval.txt > data1.txt
      fi
  ((i++))
  done
done

paste --delimiters=, data1.txt data2.txt > data.txt

rm data1.txt
rm data2.txt

