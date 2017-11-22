#!/usr/bin/env bash
dir=/home/danielbord/CLionProjects/AutomaticSolver/pipeline/points_correspondence/B
find $dir -name "points_*_left" | sed 's/^\(.*\)$/--input1 \1/g' | sort > list.left
find $dir -name "points_*_right" | sed 's/^\(.*\)$/--input2 \1/g' | sort > list.right
find $dir -name "points_*_left" | grep -o '[0-9]\+' | sed 's/^\(.*\)$/--inliers_f \1 --f_f \1.f --f_l \1.l/g' | sort > numbers.txt
paste list.left list.right numbers.txt | sed 's#^\(.*\)$#./pipeline/AutomaticSolver \1#g' | parallel

