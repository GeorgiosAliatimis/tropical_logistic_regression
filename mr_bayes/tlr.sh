#!/bin/bash

folder=$1

for dir in $folder/*/
do
	echo $dir
	cd $dir
	Rscript ../../trees_to_vectors.R 
	cd ../..
	Rscript test_case.R $dir primates_iterations.txt 0.3
	# cd $dir
	# rm v1.csv v2.csv
	# cd ../..
done
