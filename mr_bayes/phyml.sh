#!/bin/bash

for i in {10..29}
do
    ./routine.sh ../phyml_30/0$i/gene_$i phyml$i 25000 20 5000 50
done
