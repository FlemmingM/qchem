#!/bin/sh
### queue name
#BSUB -q hpc
### job name
#BSUB -J vqe_performance
### number of cores
#BSUB -n 24
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- specify that we need 4GB of memory per core/slot -- 
#BSUB -R "rusage[mem=8GB]"
### -- set walltime limit: hh:mm -- 
#BSUB -W 48:00
### Output file and error file
#BSUB -o vqe_performance.out
#BSUB -e vqe_performance.err


# Get used CPU
lscpu

# Run each molecule: minimal and Full CI

## Benzene
~/micromamba/envs/qrunch/bin/python compare_vqe_methods_fast.py benzene 6 6 40 -b STO-3G -o compare_vqe_methods
