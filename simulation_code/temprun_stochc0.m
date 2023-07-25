%% Check equivalency condition for high c0 at boundary
clear;clc
load('automated_runs/parameters/default_params_struct.mat')

params.deltaE = [0;0.05];
params.log10c0 = 2;
params.tf = params.tf;
params.log10c0_var = 1;
params.rho_at_t0 = [0.5;0.5];
%params.Omega = [0;params.deltaE(2)/params.tf*log(1+(10.^params.log10c0)/params.rho0)];
params.Omega = [0;(params.deltaE(2)./params.E(1)).*log(1+(10.^params.log10c0))/(params.rho0.*params.tf)];
params.environment_type = 'stochastic-c0';
%params.environment_type = 'deterministic';
params.max_batches = 1000;
output = serialdil_odesolver(params, 1);
%0.032+
