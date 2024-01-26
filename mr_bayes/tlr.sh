#!/bin/bash

folder=$1
iter_file=$2

for dir in $folder/*/
do
	echo $dir
	Rscript mr_bayes/trees_to_vectors.R $dir
	Rscript mr_bayes/test_case.R $dir $iter_file 0.3
	rm $dir/v1.csv $dir/v2.csv
done
