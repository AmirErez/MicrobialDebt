function [output,params] = run_random_2nutrient_community(m,log10c0,rho0,debt_type,environment_type)

load('automated_runs/parameters/default_params_struct.mat')

params.m = m;
params.log10c0 = log10c0;
params.rho0 = rho0;
params.rho_at_t0 = zeros(params.m,1) + params.rho0/params.m;
params.E = ones(params.m,1);
params.deltaE = 2*rand(params.m,1) - 1;

% Two nutrient
% params.p = 2;
% params.alpha = rand(params.m,1);
% params.alpha = sort(params.alpha,'ascend');
% params.alpha = [params.alpha, 1-params.alpha];
% params.P = rand(1);
% params.P = [params.P;1-params.P];

% One nutrient
params.p = 1;
params.alpha = ones(params.m,1);
params.P = 1;

params.environment_type = environment_type;
params.max_batches = 100000;
params.K = lognrnd(1,1,params.m,params.p);

if strcmp(debt_type,'continuous')
    params.omega = (params.E+params.deltaE).*rand(params.m,1);
    params.Omega = zeros(params.m,1);
elseif strcmp(debt_type,'punctuated')
    params.omega = zeros(params.m,1);
    params.Omega = zeros(params.m,1);
    %params.Omega = rand(params.m,1);
elseif strcmp(debt_type,'both')
    params.omega = (params.E+params.deltaE).*params.E.*rand(params.m,1);
    params.Omega = rand(params.m,1);
else
    disp('Invalid debt type.')
end

plt = 0;
output = serialdil_odesolver(params,plt);    

end

