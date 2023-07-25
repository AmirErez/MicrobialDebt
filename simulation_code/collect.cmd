#!/bin/bash
#SBATCH -N 1 # node count
#SBATCH -c 1
#SBATCH -t 12:00:00
#SBATCH --mem=4000
##SBATCH --ntasks-per-core=2
#SBATCH -o logs/collect_%A.out

module load matlab/2021a mcc
export MCR_CACHE_ROOT=/tmp/$SLURM_JOB_ID

# Stochastic runs - with timecourses
#srun matlab -novjm -r "collect_stochastic_runs('hpc_high_Omega_coexistence_10K', 4999,2,[1,4e5],100); quit;"
#srun matlab -novjm -r "collect_stochastic_runs('hpc_high_Omega_coexistence_10K_part2', 4999,2,[1,4e5],100); quit;"
#srun matlab -novjm -r "collect_precollected_stochastic_runs('hpc_high_Omega_coexistence_10K_part2', 4999,2,[1,4e5],100); quit;"
#srun matlab -novjm -r "collect_stochastic_runs('hpc_bernoulli_stochastic_coexistence',999,1,[1,2e5],100); quit;"
#srun matlab -novjm -r "collect_stochastic_runs('hpc_deterministic_stochastic_IC_sweep',1200,1,[1,3e5],100); quit;"
#srun matlab -novjm -r "collect_stochastic_runs('hpc_moderate_Omega_stochastic',999,1,[1,3e5],100); quit;"
#srun matlab -novjm -r "collect_precollected_stochastic_runs('hpc_high_Omega_fifth_noise_10K',999,1,[1,4e5],100); quit;"
#srun matlab -novjm -r "collect_precollected_stochastic_runs('hpc_tf_noise_deviation_check', 999, 1, [1,3e5], 100); quit;"
srun matlab -novjm -r "collect_precollected_stochastic_runs('hpc_tf_noise_different_starting', 2999, 1, [1,3e5], 100); quit;"

# Without timecourses
#srun matlab -novjm -r "collect_automated_serialdil('automated_runs/parameters/production_variable_K_sweep.csv', 'hpc_production_variable_K_sweep',82082, 1); quit;"
