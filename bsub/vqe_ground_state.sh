#!/bin/sh
### queue name
#BSUB -q hpc
### job name
#BSUB -J vqe_ground_state
### number of cores
#BSUB -n 1
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- specify that we need 4GB of memory per core/slot -- 
#BSUB -R "rusage[mem=10GB]"
### -- set walltime limit: hh:mm -- 
#BSUB -W 72:00
### Output file and error file
#BSUB -o vqe_ground_state.out
#BSUB -e vqe_ground_state.err


# Get used CPU
lscpu


# LiH - STO-3G, Active Space (2,2)
~/micromamba/envs/qrunch/bin/python compare_ground_state.py LiH 2 2 -m fast_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py LiH 2 2 -m oofast_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py LiH 2 2 -m adapt_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py LiH 2 2 -m ooadapt_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py LiH 2 2 -m fast_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py LiH 2 2 -m oofast_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py LiH 2 2 -m adapt_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py LiH 2 2 -m ooadapt_vqe -b STO-3G -o vqe_ground_state

# LiH - cc-pVDZ, Active Space (2,2) (Optional expansion)
~/micromamba/envs/qrunch/bin/python compare_ground_state.py LiH 2 2 -m fast_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py LiH 2 2 -m oofast_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py LiH 2 2 -m adapt_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py LiH 2 2 -m ooadapt_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py LiH 2 2 -m fast_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py LiH 2 2 -m oofast_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py LiH 2 2 -m adapt_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py LiH 2 2 -m ooadapt_vqe -b cc-pvdz -o vqe_ground_state


# H20 - STO-3G (Minimal test)
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 2 2 -m fast_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 2 2 -m oofast_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 2 2 -m adapt_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 2 2 -m ooadapt_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 2 2 -m fast_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 2 2 -m oofast_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 2 2 -m adapt_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 2 2 -m ooadapt_vqe -b STO-3G -o vqe_ground_state

~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 4 4 -m fast_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 4 4 -m oofast_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 4 4 -m adapt_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 4 4 -m ooadapt_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 4 4 -m fast_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 4 4 -m oofast_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 4 4 -m adapt_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 4 4 -m ooadapt_vqe -b STO-3G -o vqe_ground_state

~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 6 5 -m fast_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 6 5 -m oofast_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 6 5 -m adapt_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 6 5 -m ooadapt_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 6 5 -m fast_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 6 5 -m oofast_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 6 5 -m adapt_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 6 5 -m ooadapt_vqe -b STO-3G -o vqe_ground_state

~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 8 6 -m fast_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 8 6 -m oofast_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 8 6 -m adapt_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 8 6 -m ooadapt_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 8 6 -m fast_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 8 6 -m oofast_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 8 6 -m adapt_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 8 6 -m ooadapt_vqe -b STO-3G -o vqe_ground_state

~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 10 7 -m fast_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 10 7 -m oofast_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 10 7 -m adapt_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 10 7 -m ooadapt_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 10 7 -m fast_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 10 7 -m oofast_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 10 7 -m adapt_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H20 10 7 -m ooadapt_vqe -b STO-3G -o vqe_ground_state

# H20 - cc-pvdz (Relevant for qLR shot noise studies)
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 4 4 -m fast_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 4 4 -m oofast_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 4 4 -m adapt_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 4 4 -m ooadapt_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 4 4 -m fast_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 4 4 -m oofast_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 4 4 -m adapt_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 4 4 -m ooadapt_vqe -b cc-pvdz -o vqe_ground_state

~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 6 5 -m fast_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 6 5 -m oofast_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 6 5 -m adapt_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 6 5 -m ooadapt_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 6 5 -m fast_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 6 5 -m oofast_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 6 5 -m adapt_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 6 5 -m ooadapt_vqe -b cc-pvdz -o vqe_ground_state

~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 8 6 -m fast_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 8 6 -m oofast_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 8 6 -m adapt_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 8 6 -m ooadapt_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 8 6 -m fast_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 8 6 -m oofast_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 8 6 -m adapt_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 8 6 -m ooadapt_vqe -b cc-pvdz -o vqe_ground_state

~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 4 4 -m fast_vqe_as -b cc-pvtz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 4 4 -m oofast_vqe_as -b cc-pvtz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 4 4 -m adapt_vqe_as -b cc-pvtz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 4 4 -m ooadapt_vqe_as -b cc-pvtz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 4 4 -m fast_vqe -b cc-pvtz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 4 4 -m oofast_vqe -b cc-pvtz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 4 4 -m adapt_vqe -b cc-pvtz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py H2O 4 4 -m ooadapt_vqe -b cc-pvtz -o vqe_ground_state

# BeH2
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 4 -m fast_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 4 -m oofast_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 4 -m adapt_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 4 -m ooadapt_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 4 -m fast_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 4 -m oofast_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 4 -m adapt_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 4 -m ooadapt_vqe -b STO-3G -o vqe_ground_state

~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 6 -m fast_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 6 -m oofast_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 6 -m adapt_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 6 -m ooadapt_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 6 -m fast_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 6 -m oofast_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 6 -m adapt_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 6 -m ooadapt_vqe -b STO-3G -o vqe_ground_state

~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 6 7 -m fast_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 6 7 -m oofast_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 6 7 -m adapt_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 6 7 -m ooadapt_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 6 7 -m fast_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 6 7 -m oofast_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 6 7 -m adapt_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 6 7 -m ooadapt_vqe -b STO-3G -o vqe_ground_state

~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 4 -m fast_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 4 -m oofast_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 4 -m adapt_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 4 -m ooadapt_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 4 -m fast_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 4 -m oofast_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 4 -m adapt_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 4 -m ooadapt_vqe -b cc-pvdz -o vqe_ground_state

~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 6 -m fast_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 6 -m oofast_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 6 -m adapt_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 6 -m ooadapt_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 6 -m fast_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 6 -m oofast_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 6 -m adapt_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 4 6 -m ooadapt_vqe -b cc-pvdz -o vqe_ground_state

~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 6 7 -m fast_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 6 7 -m oofast_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 6 7 -m adapt_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 6 7 -m ooadapt_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 6 7 -m fast_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 6 7 -m oofast_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 6 7 -m adapt_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py BeH2 6 7 -m ooadapt_vqe -b cc-pvdz -o vqe_ground_state


# Benzene
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 6 6 -m fast_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 6 6 -m oofast_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 6 6 -m adapt_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 6 6 -m ooadapt_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 6 6 -m fast_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 6 6 -m oofast_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 6 6 -m adapt_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 6 6 -m ooadapt_vqe -b STO-3G -o vqe_ground_state

~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 8 8 -m fast_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 8 8 -m oofast_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 8 8 -m adapt_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 8 8 -m ooadapt_vqe_as -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 8 8 -m fast_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 8 8 -m oofast_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 8 8 -m adapt_vqe -b STO-3G -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 8 8 -m ooadapt_vqe -b STO-3G -o vqe_ground_state

~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 6 6 -m fast_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 6 6 -m oofast_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 6 6 -m adapt_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 6 6 -m ooadapt_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 6 6 -m fast_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 6 6 -m oofast_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 6 6 -m adapt_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 6 6 -m ooadapt_vqe -b cc-pvdz -o vqe_ground_state

~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 8 8 -m fast_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 8 8 -m oofast_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 8 8 -m adapt_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 8 8 -m ooadapt_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 8 8 -m fast_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 8 8 -m oofast_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 8 8 -m adapt_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 8 8 -m ooadapt_vqe -b cc-pvdz -o vqe_ground_state

~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 6 6 -m fast_vqe_as -b cc-pvtz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 6 6 -m oofast_vqe_as -b cc-pvtz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 6 6 -m adapt_vqe_as -b cc-pvtz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 6 6 -m ooadapt_vqe_as -b cc-pvtz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 6 6 -m fast_vqe -b cc-pvtz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 6 6 -m oofast_vqe -b cc-pvtz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 6 6 -m adapt_vqe -b cc-pvtz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 6 6 -m ooadapt_vqe -b cc-pvtz -o vqe_ground_state

~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 8 8 -m fast_vqe_as -b cc-pvtz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 8 8 -m oofast_vqe_as -b cc-pvtz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 8 8 -m adapt_vqe_as -b cc-pvtz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 8 8 -m ooadapt_vqe_as -b cc-pvtz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 8 8 -m fast_vqe -b cc-pvtz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 8 8 -m oofast_vqe -b cc-pvtz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 8 8 -m adapt_vqe -b cc-pvtz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py benzene 8 8 -m ooadapt_vqe -b cc-pvtz -o vqe_ground_state

# R-methyloxirane
~/micromamba/envs/qrunch/bin/python compare_ground_state.py R-methyloxirane 6 6 -m fast_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py R-methyloxirane 6 6 -m oofast_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py R-methyloxirane 6 6 -m adapt_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py R-methyloxirane 6 6 -m ooadapt_vqe_as -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py R-methyloxirane 6 6 -m fast_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py R-methyloxirane 6 6 -m oofast_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py R-methyloxirane 6 6 -m adapt_vqe -b cc-pvdz -o vqe_ground_state
~/micromamba/envs/qrunch/bin/python compare_ground_state.py R-methyloxirane 6 6 -m ooadapt_vqe -b cc-pvdz -o vqe_ground_state
