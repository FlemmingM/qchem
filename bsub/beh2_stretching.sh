#!/bin/sh
### queue name
#BSUB -q hpc
### job name
#BSUB -J beh2_stretching
### number of cores
#BSUB -n 1
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- specify that we need 4GB of memory per core/slot -- 
#BSUB -R "rusage[mem=10GB]"
### -- set walltime limit: hh:mm -- 
#BSUB -W 48:00
### Output file and error file
#BSUB -o beh2_stretching.out
#BSUB -e beh2_stretching.err


# Get used CPU
lscpu
~/micromamba/envs/qrunch/bin/python beh2_stretching.py
