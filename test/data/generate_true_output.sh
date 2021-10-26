#!/bin/bash

conda activate decifer_dev
decifer ./mutations.tsv -p ./purity.tsv -k 5 -K 8 -r 20 --seed 17 -j 4  
