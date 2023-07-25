function [tf,Qf] = compute_single_species_tf(Q,log10c0,rho0,E,gtype,h)

%This function computes tf given a Q

params.p = 1;
params.m = 1;
params.K = 1;
params.batch_length = Inf;
params.P = 1;
params.Omega = 0;
params.omega = 0;
params.gtype = gtype;
params.h = h;
params.deltaE = 0;
params.alpha = 1;
params.rho0 = rho0;
params.log10c0 = log10c0;
params.E = E;
params.rho_at_t0 = 1;
params.Q = Q;


y0 = [10.^params.log10c0; params.rho_at_t0];
tspan=[0,Inf];
options = odeset('NonNegative',1,'RelTol',1e-11, ...
    'Events', @(t,y)eventfun_Q(t,y,params,Q));%, ...

[t,y] = ode15s(@(t,y) odefun(t,y,params), tspan, y0, options);

tf = t(end);
Qf = (10^params.log10c0 - y(end,1))./(10^params.log10c0);

end