function dydt = odefun(t,y,params)
% In the form of dy/dt = f(y). Include both c_i and rho_sigma.

alpha_eff = params.alpha.*(repmat(params.E + params.deltaE,1,params.p));

if strcmp(params.gtype,'hill')
    gfun = @(x,K) ((x.^params.h)./((K.^params.h) + (x.^params.h)));
elseif strcmp(params.gtype,'linear')
    gfun = @(x,K) x./K;
else
    disp('Invalid growth function type')
end

c_i = y(1:params.p);
rho_sigma=y((params.p+1):end);

alpha_g_rho_mat = alpha_eff.*gfun(repmat(c_i',params.m,1),params.K).*repmat(rho_sigma,1,params.p);

dydt = zeros(length(c_i)+length(rho_sigma),1);
dydt(1:params.p) = -sum(alpha_g_rho_mat,1)';
dydt((params.p+1):end) = sum(alpha_g_rho_mat,2)  - params.omega.*rho_sigma;

end

