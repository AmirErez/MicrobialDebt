# MicrobialDebt
Code and data for the manuscript "Grow now, pay later: when should a bacterium go into debt?"

* simulation_code/ contains the code for the deterministic and stochastic 
simulations.
* simulation_code/automated_runs/results/ contains data used for figure generation. The md5.txt checksum can be used to verify data integrity.
* Note that simulation_code/automated_runs/results/hpc_high_Omega_coexistence_10K_part1_collected_results.mat.* needs to be merged to a single file,
it was split to overcome github file size limits.
* simulation_code/fig_*.m are the scripts for plotting the figures in 
the manuscript.
* data_fitting/ fits the chemostat dynamics from Notley-McRobb et al.
* figures/ contains the figures used in the manuscript, produced by 
simulation_code/fig_*.m

