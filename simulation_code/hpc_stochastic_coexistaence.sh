#!/bin/bash
##SBATCH --array=3000-3000
##SBATCH --array=1-1000
#SBATCH --array=1-1
##SBATCH --array=1-5000
#SBATCH -N 1 # node count
#SBATCH --overcommit
#SBATCH --ntasks-per-core=2
#SBATCH -c 1
#SBATCH -t 23:00:00
#SBATCH --mem=4000
#SBATCH -o logs/K_coexistence_%A_%a.out

# sent 37000
actual=$(($SLURM_ARRAY_TASK_ID+9999))
#actual=$(($SLURM_ARRAY_TASK_ID))



# table_file='automated_runs/parameters/check_stochastic_coexistence.csv'
# results_folder='automated_runs/results/hpc_check_stochastic_coexistence/'
# nrun=1

# table_file='automated_runs/parameters/critical_stochastic_coexistence.csv'
# results_folder='automated_runs/results/hpc_critical_stochastic_coexistence/'
# nrun=1

# table_file='automated_runs/parameters/opposite_stochastic_coexistence.csv'
# results_folder='automated_runs/results/hpc_opposite_stochastic_coexistence/'
# nrun=2

#table_file='automated_runs/parameters/high_Omega_coexistence_10K.csv'
#results_folder='automated_runs/results/hpc_high_Omega_coexistence_10K_part2/'
#nrun=2

# table_file='automated_runs/parameters/bernoulli_stochastic_coexistence.csv'
# results_folder='automated_runs/results/hpc_bernoulli_stochastic_coexistence/'
# nrun=1


#table_file='automated_runs/parameters/deterministic_stochastic_IC_sweep.csv'
#results_folder='automated_runs/results/hpc_deterministic_stochastic_IC_sweep/'
#nrun=1

#table_file='automated_runs/parameters/moderate_Omega_stochastic.csv'
#results_folder='automated_runs/results/hpc_moderate_Omega_stochastic/'
#nrun=1

#table_file='automated_runs/parameters/high_Omega_fifth_noise_10K.csv'
#results_folder='automated_runs/results/hpc_high_Omega_fifth_noise_10K'
#nrun=1

table_file='automated_runs/parameters/production_variable_K_sweep.csv'
results_folder='automated_runs/results/hpc_production_variable_K_sweep'
nrun=1

#table_file='automated_runs/parameters/tf_noise_deviation_check.csv'  # 3000 runs
#results_folder='automated_runs/results/hpc_tf_noise_deviation_check'
#nrun=1

#table_file='automated_runs/parameters/tf_noise_different_starting.csv' # 3000 runs
#results_folder='automated_runs/results/hpc_tf_noise_different_starting'
#nrun=1


mkdir -p $results_folder
default_params_file='automated_runs/parameters/default_params_struct.mat'
results_file="$results_folder/sim_group_${nrun}.mat"
if [ -e $results_file ]; then
    echo "File $results_file exists, not overwriting, exiting."
    exit 0
fi
str="run_automated_run_serialdil.sh /usr/local/matlab/default $default_params_file $table_file $actual $results_folder $nrun"
echo $str
$str
#matlab -nojvm -r "addpath(genpath('/home/jaime/microbial_debt/')); automated_run_serialdil( '$table_file' , $SLURM_ARRAY_TASK_ID , '$results_folder',1); quit();"

