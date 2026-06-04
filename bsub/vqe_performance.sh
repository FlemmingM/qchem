#!/bin/sh
### queue name
#BSUB -q hpc
### job name
#BSUB -J vqe_performance
### number of cores
#BSUB -n 1
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- specify that we need 4GB of memory per core/slot -- 
#BSUB -R "rusage[mem=10GB]"
### -- set walltime limit: hh:mm -- 
#BSUB -W 72:00
### Output file and error file
#BSUB -o vqe_performance.out
#BSUB -e vqe_performance.err


# Get used CPU
lscpu

# Run each molecule: minimal and Full CI

## LiH - Done
# ~/micromamba/envs/qrunch/bin/python compare_vqe_methods.py LiH 2 2 3 -b STO-3G -o compare_vqe_methods
# ~/micromamba/envs/qrunch/bin/python compare_vqe_methods.py LiH 4 6 40 -b STO-3G -o compare_vqe_methods

## H2O
~/micromamba/envs/qrunch/bin/python compare_vqe_methods.py H2O 2 2 3 -b STO-3G -o compare_vqe_methods
~/micromamba/envs/qrunch/bin/python compare_vqe_methods.py H2O 4 4 20 -b STO-3G -o compare_vqe_methods
~/micromamba/envs/qrunch/bin/python compare_vqe_methods.py H2O 6 5 40 -b STO-3G -o compare_vqe_methods
~/micromamba/envs/qrunch/bin/python compare_vqe_methods.py H2O 8 6 40 -b STO-3G -o compare_vqe_methods
~/micromamba/envs/qrunch/bin/python compare_vqe_methods.py H2O 10 7 40 -b STO-3G -o compare_vqe_methods

## BeH2
~/micromamba/envs/qrunch/bin/python compare_vqe_methods.py BeH2 4 4 20 -b STO-3G -o compare_vqe_methods
~/micromamba/envs/qrunch/bin/python compare_vqe_methods.py BeH2 6 7 40 -b STO-3G -o compare_vqe_methods

## R-methyloxirane
~/micromamba/envs/qrunch/bin/python compare_vqe_methods.py R-methyloxirane 6 6 40 -b STO-3G -o compare_vqe_methods_ecd
