%This script generates the automated parameter struct used in the automated
%run code

clear;clc

params.p = 1;
params.m = 2;
params.K = 1;
params.max_batches = 1e6;
params.tf = compute_single_species_tf(0.99,4,1,1);
params.E = [1;1];
params.P = 1;
params.rho0 = 1;
params.rho_at_t0 = [0.5*params.rho0;0.5*params.rho0];
params.deltaE = [0;0];
params.gtype = 'hill';
params.h = 1;
params.omega = [0;0];
params.Omega = [0;0];
params.log10c0 = 0;
params.alpha = [1;1];
params.errtype = 1;
params.environment_type = 'deterministic';
params.log10c0_var = 1;
params.tf_var = 1;

save('automated_runs/parameters/default_params_struct.mat')
