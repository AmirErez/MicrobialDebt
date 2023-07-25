K=1;
rho0=1;
log10c0=2;


% Calculate and plot
params.p = 1;
params.m = 1;
params.K = K;
params.tf = compute_single_species_tf(0.99,4,1,1,'hill',1);
params.E = 1;
params.P = 1;
params.rho0 = rho0;
params.rho_at_t0 = params.rho0;
params.deltaE = 0;
params.Omega = 0;
params.omega = 0;
params.alpha = 1;
params.errtype = 1;
params.log10c0 = log10c0;
params.dilution_method = 'floating';
params.gtype = 'hill';
params.h = 1;
params.environment_type = 'deterministic';
plt = 0;

gfun_sq = @(x,K) x./(x + K).^2;

deltaEs = 0:0.05:1;

newfigure(3,2);
% figure(666);
hold on
omap = colormap(copper(length(deltaEs)));
for ii = 1:length(deltaEs)
    c0=10.^params.log10c0;
    params.deltaE = deltaEs(ii);
    %params.omega = 1/params.tf*(log(params.rho0/(params.rho0+c0)) + log(1+c0/params.rho0*(1+params.deltaE/params.E)));
    params.omega = params.deltaE/params.tf/params.E*log(1+c0/params.rho0);
    [rho_sigma, Nr, c_i, t] = multispeciesbatch_odesolver(params,0);
    plot(t/params.tf,rho_sigma,'-','Color',omap(ii,:), 'LineWidth',1.5)
end
ylim([0.5,100])
c=colorbar();
c.Ticks=[0,1];
% c.TickLabels={'$\Delta_E=0,\\ \omega=0$', '$\Delta_E=1,\\ \omega=0.50$'};
c.TickLabels={'$\begin{array}{l} \Delta_E=0 \\ \omega=0 \end{array}$', '$\begin{array}{c} \Delta_E=1 \\ \omega=0.5 \end{array}$'}
c.TickLabelInterpreter='Latex';
% c.Location = 'southeast';
set(gca,'YScale','log');
xlabel('$t/t_f$', 'Interpreter','latex');
ylabel('Single species, $\rho/\rho_0$', 'Interpreter','latex');

%% Now plot the two species coexisting curve



params.m = 2;
params.E = [1;1];
params.rho_at_t0 = [0.5*params.rho0;0.5*params.rho0];
params.deltaE = [1;0];
params.omega = [0.33;0];
params.Omega = [0;0];
params.alpha = [1;1];
[rho_sigma, Nr, c_i, t] = multispeciesbatch_odesolver(params,0);

ax=axes();
ax.Position = [0.3194    0.2569    0.2686    0.2778];
hold on
debtor_color = [206 37 123]/255;
nondebtor_color = [15 104 194]/255;
plot(t/params.tf,rho_sigma(1,:),'LineWidth', 1.5, 'Color',debtor_color); 
plot(t/params.tf,rho_sigma(2,:),'LineWidth', 1.5,'Color',nondebtor_color); 
set(gca,'YScale','log')
xticks([])
yticks([])
ylim([0.5,100])
box("on")


print(gcf,'-dpng', '../figures/fig_single_species_debtor.png', '-r600')
print(gcf,'-dsvg', '../figures/fig_single_species_debtor.svg')