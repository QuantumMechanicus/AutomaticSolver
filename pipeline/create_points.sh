#!/usr/bin/env bash
dir=/home/danielbord/CLionProjects/AutomaticSolver/pipeline/points_correspondence/B
out_dir=./automatic_solver_results_test/
list_l=${out_dir}list.left
list_r=${out_dir}list.right
numbers=${out_dir}numbers.txt
mkdir -p $out_dir
find $dir -name "points_*_left" | sed 's/^\(.*\)$/--input1 \1/g' | sort > $list_l
find $dir -name "points_*_right" | sed 's/^\(.*\)$/--input2 \1/g' | sort > $list_r
find $dir -name "points_*_left" | grep -o '[0-9]\+' | sed 's/^\(.*\)$/--inliers_f \1 --fund_f \1.f --lambd_f \1.l/g' | sort > $numbers
echo $numbers
#paste $list_l $list_r $numbers | sed 's#^\(.*\)$#./GroebnerAutomaticSolver \1#g' | parallel

