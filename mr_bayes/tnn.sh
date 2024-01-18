#!/bin/bash

folder=$1

cd $folder
for i in {0..99}
do
    cd $i
    python3 ~/Downloads/primates_convergence/tnn.py 50
    cd ..
done
