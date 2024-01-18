#!/bin/bash

# Arguments 
# $1: gene file (nexus) studied
# $2: name of directory to be created (or extended)
# $3: ngen; make sure that this is an integer multiple of diagnostic frequency
# $4: sampling frequency
# $5: diagnostic frequency; make sure that this is an integer multiple of sampling frequency
# $6: the number of times AUCs and ASDSFs are competed, using independent MCMC runs  

EXTEND="true"

if [ $EXTEND == "false" ] 
then 
    cp $1 gene
    rm -rf $2
    mkdir $2
    mv gene $2
    cd $2
    python3 ../create_mb_file.py gene $3 $4 $5
elif [ $EXTEND == "true" ]
then
    cd $2 
fi

# start=`find . -maxdepth 1 -type d | wc -l`
start=1
start=$(( start - 1 ))
end=$(( start + $6 -1))

for (( i=$start; i<=$end; i++ ))
do
    printf "no" | mb -i mb_file > out
    grep "Average standard deviation of split frequencies" out |  grep -Eo '[0-9\.]+'>> asdsfs

    #mkdir trees1 trees2
    mv gene.run1.t t1
    mv gene.run2.t t2
    #Rscript ../trees_to_vectors.R $4 $5
    rm gene.* out
    
    #python3 ../tnn.py $7
    
    #mkdir $i
    
    mv t1 $i
    mv t2 $i
    mv asdsfs $i
    #mv aucs $i
    #mv trees1 $i
    #mv trees2 $i
done
cd ..