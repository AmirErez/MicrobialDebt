function [Nr_den_sq, Nrs] = compute_Nr_den_sq(log10c0_range,rho0,K,deltaE)


params.p = 1;
params.m = 1;
params.K = K;
params.tf = compute_single_species_tf(0.99,4,1,1,'hill',1);
params.E = 1;
params.P = 1;
params.rho0 = rho0;
params.rho_at_t0 = params.rho0;
params.deltaE = deltaE;
params.Omega = 0;
params.omega = 0;
params.alpha = 1;
params.errtype = 1;
params.dilution_method = 'floating';
params.gtype = 'hill';
params.h = 1;
params.environment_type = 'deterministic';
plt = 0;

gfun_sq = @(x,K) x./((x + K).^2);

for i = 1:length(log10c0_range)
    
    params.log10c0 = log10c0_range(i);
    [rho_sigma, Nr, c_i, t] = multispeciesbatch_odesolver(params,plt);
%     Nr_den_sq(i) = (params.E  + params.deltaE).*trapz(t,transpose(gfun_sq(c_i,K))); %Compute integral of growth
    Nr_den_sq(i) = trapz(t,transpose(gfun_sq(c_i,K))); %Compute integral of growth
    Nrs(i) = Nr;
end


end

