#!/usr/bin/bash

# This code uses the split-apply-combine template. 
# The main matlab script is compiled to avoid license problems when running >10 jobs at once
# Amir Erez 2022-06-17
# amir.erez1@mail.huji.ac.il

# Compile
sbatch compile.cmd
# Run
sbatch hpc_check_stochastic_coexistaence.sh

#Collect
srun matlab -novjm -r "collect_stochastic_runs('hpc_critical_stochastic_coexistence', 999,1,[1,2e5],100); quit;"
