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

## H2O
# ~/micromamba/envs/qrunch/bin/python compare_vqe_methods.py H2O 4 4 20 -b cc-pvdz -o compare_vqe_methods

~/micromamba/envs/qrunch/bin/python compare_vqe_methods.py H2O 6 5 40 -b cc-pvdz -o compare_vqe_methods
~/micromamba/envs/qrunch/bin/python compare_vqe_methods.py H2O 8 6 40 -b cc-pvdz -o compare_vqe_methods
