#!/bin/bash

Rscript --vanilla density.r Distr
est_lambda=`cat estimated_lambda`
python undistort.py --coeff "$est_lambda" --img 'D.JPG' --output 'D'