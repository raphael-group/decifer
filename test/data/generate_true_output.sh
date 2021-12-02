#!/bin/bash

conda activate decifer_dev
decifer \
/Users/bjarnold/Princeton_DataX/decifer/test/data2/decifer_input.tsv \
--purityfile \
/Users/bjarnold/Princeton_DataX/decifer/test/data2/purity.tsv \
--betabinomial \
--segfile \
/Users/bjarnold/Princeton_DataX/decifer/test/data2/segfile.tsv \
--snpfile \
/Users/bjarnold/Princeton_DataX/decifer/test/data2/snpfile.tsv \
--seed \
17 \
--mink \
6 \
--maxk \
6 \
--debug \
