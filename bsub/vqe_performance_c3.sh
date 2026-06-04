#!/bin/sh
### queue name
#BSUB -q hpc
### job name
#BSUB -J vqe_performance
### number of cores
#BSUB -n 3
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- specify that we need 4GB of memory per core/slot -- 
#BSUB -R "rusage[mem=10GB]"
### -- set walltime limit: hh:mm -- 
#BSUB -W 24:00
### Output file and error file
#BSUB -o compare_vqe_methods_c3.out
#BSUB -e compare_vqe_methods_c3.err


# Get used CPU
lscpu

## H2O
~/micromamba/envs/qrunch/bin/python compare_vqe_methods_simple.py H2O 4 4 20 -b STO-3G -o compare_vqe_methods_c3
