#!/usr/bin/env bash
dir=/home/danielbord/CLionProjects/AutomaticSolver/pipeline/points_correspondence/T
out_dir=../automatic_solver_results_T/
list_l=${out_dir}list.left
list_r=${out_dir}list.right
find $dir -name "points_*_left" | sed 's/^\(.*\)$/\1/g' | sort > $list_l
find $dir -name "points_*_right" | sed 's/^\(.*\)$/\1/g' | sort > $list_r
rm -f ${out_dir}all.ll
rm -f ${out_dir}all.ff
find ./$out_dir -name "*.l" | sort | xargs -n 1 -I^ cat ^ >> ${out_dir}all.ll
find ./$out_dir -name "*.f" | sort | xargs -n 1 -I^ cat ^ >> ${out_dir}all.ff
ones=`find $dir -name "points_*_left" | sort | grep -o '[0-9]\+' | sed 's/^\(.*\)$/\1/g' | xargs -n 1 printf ${out_dir}'%s_left '`
twos=`find $dir -name "points_*_left" | sort | grep -o '[0-9]\+' | sed 's/^\(.*\)$/\1/g' | xargs -n 1 printf ${out_dir}'%s_right '`

s_left=$(paste -s $list_l)
s_right=$(paste -s $list_r)
s_left="--left_corr_f "${s_left}
s_right="--right_corr_f "${s_right}
../NonLinearOptimizator --n_iters 3 $s_left $s_right --n_pic 4 --w 5184 --h 3456 --fund_f ${out_dir}all.ff --lambda_f ${out_dir}all.ll --nlambda 2 --q 0.1
