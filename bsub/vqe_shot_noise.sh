#!/bin/sh
### queue name
#BSUB -q hpc
### job name
#BSUB -J vqe_shot_noise
### number of cores
#BSUB -n 24
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- specify that we need 4GB of memory per core/slot -- 
#BSUB -R "rusage[mem=8GB]"
### -- set walltime limit: hh:mm -- 
#BSUB -W 48:00
### Output file and error file
#BSUB -o vqe_shot_noise.out
#BSUB -e vqe_shot_noise.err


# Get used CPU
lscpu


## H2O
~/micromamba/envs/qrunch/bin/python compare_vqe_shot_noise.py H2O 4 4 20 -b STO-3G -r 10 -o vqe_shot_noise

