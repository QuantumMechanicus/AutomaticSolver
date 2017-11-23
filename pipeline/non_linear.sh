#!/usr/bin/env bash
dir=/home/danielbord/CLionProjects/AutomaticSolver/pipeline/points_correspondence/B
out_dir=./automatic_solver_results/
list_l=${out_dir}list.left
list_r=${out_dir}list.right
find $dir -name "points_*_left" | sed 's/^\(.*\)$/--input1 \1/g' | sort > $list_l
find $dir -name "points_*_right" | sed 's/^\(.*\)$/--input2 \1/g' | sort > $list_l
rm -f ${out_dir}all.ll
rm -f ${out_dir}all.ff
find ./ -name "*.l" | sort | xargs -n 1 -I^ cat ^ >> ${out_dir}all.ll
find ./ -name "*.f" | sort | xargs -n 1 -I^ cat ^ >> ${out_dir}all.ff
ones=`find $dir -name "points_*_left" | sort | grep -o '[0-9]\+' | sed 's/^\(.*\)$/\1/g' | xargs -n 1 printf ${out_dir}'%s_left '`
twos=`find $dir -name "points_*_left" | sort | grep -o '[0-9]\+' | sed 's/^\(.*\)$/\1/g' | xargs -n 1 printf ${out_dir}'%s_right '`
./NonLinearOptimizator --left_inl_f $ones --n_pic 9 --right_inl_f $twos --fund_f ${out_dir}all.ff --lambda_f ${out_dir}all.ll --nlambda 2
