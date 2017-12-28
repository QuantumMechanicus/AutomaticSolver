#!/usr/bin/env bash
dir=/home/danielbord/CLionProjects/AutomaticSolver/pipeline/points_correspondence/T
out_dir=../automatic_solver_results_T/
list_l=${out_dir}list.left
list_r=${out_dir}list.right
numbers=${out_dir}numbers.txt
mkdir -p $out_dir
find $dir -name "points_*_left" | sed 's/^\(.*\)$/--input1 \1/g' | sort > $list_l
find $dir -name "points_*_right" | sed 's/^\(.*\)$/--input2 \1/g' | sort > $list_r
find $dir -name "points_*_left" | grep -o '[0-9]\+' | sed "sX^\(.*\)\$X --iters 500000 --w 5184 --h 3456 --inliers_f $out_dir\1 --fund_f $out_dir\1.f --q 0.2 --lambd_f $out_dir\1.lXg" | sort > $numbers
paste $list_l $list_r $numbers | sed 's#^\(.*\)$#../GroebnerAutomaticSolver \1#g' | parallel

