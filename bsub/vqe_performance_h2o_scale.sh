#!/bin/sh
### queue name
#BSUB -q hpc
### job name
#BSUB -J vqe_performance_h2o_scale
### number of cores
#BSUB -n 1
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- specify that we need 4GB of memory per core/slot -- 
#BSUB -R "rusage[mem=10GB]"
### -- set walltime limit: hh:mm -- 
#BSUB -W 72:00
### Output file and error file
#BSUB -o vqe_performance_h2o_scale.out
#BSUB -e vqe_performance_h2o_scale.err


# Get used CPU
lscpu

~/micromamba/envs/qrunch/bin/python compare_vqe_methods.py H2O 10 7 40 -b cc-pvtz -o vqe_performance_h2o_scale
