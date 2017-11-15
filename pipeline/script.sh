#!/bin/bash
input1=$1
input2=$2
distribution=$3
image=$4
iters=${5:-10000}
threshold=${6:-5}
threshold2=${7:--5}
w=${8:-7360.0}
h=${9:-4912.0}
echo "First image detected points: $1"
echo "Second image detected points: $2"
echo "output distribution: $3"
echo "image to undistort: $4"
echo "number of iterations: $iters"
echo "Threshold: $threshold"
echo "Threshold: $threshold2"
echo "width: $w"
echo "height: $h"
echo "Run Automatic solver to generate distribution"
./AutomaticSolver --input1 "$input1" --input2 "$input2" --f "$distribution" --iters "$iters" --w "$w" --h "$h" --threshold1 "$threshold" --threshold2 "$threshold2"
echo "Run Automatic density estimation"
Rscript --vanilla density.r "$distribution"
echo "Run undistortion"
est_lambda=`cat estimated_lambda`
python undistort.py --coeff "$est_lambda" --img "$image"