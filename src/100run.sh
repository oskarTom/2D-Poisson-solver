#!/bin/bash

i=1
sum=0
while [ $i -le 10 ]
do
  time=$(mpirun -np 4 ../run/run)
  echo "${time}"
  i=$(($i+1))
  sum=$(awk "BEGIN {print $sum+$time; exit}")
done
sum=$(awk "BEGIN {print $sum/10; exit}")
echo $sum
