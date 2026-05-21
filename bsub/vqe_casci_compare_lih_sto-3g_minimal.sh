#!/bin/sh
### queue name
#BSUB -q hpc
### job name
#BSUB -J vqe_casci_compare
### number of cores
#BSUB -n 4
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- specify that we need 4GB of memory per core/slot -- 
#BSUB -R "rusage[mem=8GB]"
### -- set walltime limit: hh:mm -- 
#BSUB -W 24:00
### Output file and error file
#BSUB -o vqe_casci_compare_%J.out
#BSUB -e vqe_casci_compare_%J.err


~/micromamba/envs/qrunch/bin/python compare_vqe_casci.py LiH 2 2 4 -b STO-3G -o "vqe_casci_compare_lih_sto-3g_minimal_${LSB_JOBID}"
