#!/bin/sh
# This file does the same thing as qsub -V flag,
# restores LD_LIBRARY_PATH so that I can use python2.7 with qsub

source /broad/software/scripts/useuse
reuse -q python-2.7

"$@"




# python2.7 /broad/compbio/maxwshen/Kellis/1-MAKETRAINTEST/foreground_background.py

# for i in {1..1000}
# do
#   echo $i
#   sleep 1s
# done