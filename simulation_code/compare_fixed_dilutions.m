clear;clc

clear;clc

params.p = 1;
params.m = 2;
params.K = 1;
params.max_batches = 1e4;
params.batch_length = Inf;
params.E = [1;1];
params.P = 1;
params.d = [0;0];
params.rho0 = 1;
params.rho_at_t0 = [1/2;1/2];
params.psi = [1;2];
params.theta = [0;0];
params.log10c0 = 1;
params.alpha = [1;1];
params.errtype = 1;
params.Omega = [1;1];
plt = 0;

output = serialdil_odesolver(params, plt);

