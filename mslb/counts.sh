#!/bin/bash

tol=$1
for i in 0 1 2 3 4; do
  awk -v err_tol=$tol -f process_predictions.awk out | grep "phase $i" | wc -l
done
