#!/usr/bin/env bash
dir=/home/danielbord/CLionProjects/AutomaticSolver/pipeline/points_correspondence/B
find $dir -name "points_*_left" | sed 's/^\(.*\)$/--input1 \1/g' | sort > list.left
find $dir -name "points_*_right" | sed 's/^\(.*\)$/--input2 \1/g' | sort > list.right
rm -f all.ll
rm -f all.ff
find ./ -name "*.l" | sort | xargs -n 1 -I^ cat ^ >> all.ll
find ./ -name "*.f" | sort | xargs -n 1 -I^ cat ^ >> all.ff
ones=`find $dir -name "points_*_left" | sort | grep -o '[0-9]\+' | sed 's/^\(.*\)$/\1/g' | xargs -n 1 printf '%s1 '`
twos=`find $dir -name "points_*_left" | sort | grep -o '[0-9]\+' | sed 's/^\(.*\)$/\1/g' | xargs -n 1 printf '%s2 '`
./AutomaticSolver --f_names $ones --f_names2 $twos --f_f all.ff --f_l all.ll --optim 1 --nlambda 2
