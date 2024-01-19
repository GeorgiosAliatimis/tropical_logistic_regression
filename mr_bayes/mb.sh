#!/bin/bash

# Arguments 
# $1: gene file (nexus) studied
# $2: name of directory to be created
# $3: ngen; number of MCMC iterations;  make sure that this is an integer multiple of diagnostic frequency
# $4: sampling frequency
# $5: diagnostic frequency; make sure that this is an integer multiple of sampling frequency
# $6: the number of times AUCs and ASDSFs are competed, using independent MCMC runs  

cp $1 gene
mkdir $2
mv gene $2
cd $2
python3 ../mr_bayes/create_mb_file.py gene $3 $4 $5

# start=`find . -maxdepth 1 -type d | wc -l`
start=1
end=$6

for (( i=$start; i<=$end; i++ ))
do
    printf "no" | mb -i mb_file > out
    grep "Average standard deviation of split frequencies" out |  grep -Eo '[0-9\.]+'>> asdsfs

    mv gene.run1.t t1.nex
    mv gene.run2.t t2.nex
    rm gene.* out
    mkdir $i
    mv t1.nex $i
    mv t2.nex $i
    mv asdsfs $i
done
cd ..
