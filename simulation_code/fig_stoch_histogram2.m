clear

s2=load('automated_runs/results/hpc_high_Omega_fifth_noise_10K_precollected_results');
batches = s2.precollected_final.batch_cell{1};
rho1_mat = s2.precollected_final.rho1_mat;
rho2_mat = s2.precollected_final.rho2_mat;


newfigure(3,2);
set(gca,'FontSize',11)
hold on

sz=500;
nsteps=8;
sims=sz*nsteps;
smap=colormap(copper(nsteps));
sz=size(rho1_mat,2)/nsteps;

%These are the theoretical mean and std of log ratio increments ln(r_(b+1))-ln(r_b) 
%The variables used in this simulation:
c0=100 ;
rho0=1 ;
E=1 ;
deltaE=0.05 ;
Omega=0.002508737448117671 ;
delta_Omega2=0.000250873744811767 ;
max_batches=400000 ;
tf=92.0176177196183 ;
theory_mu=(Omega*tf-deltaE/E*log(1+c0/rho0));
theory_sigma=sqrt(tf^2*(delta_Omega2^2)/3);

for ii=1:nsteps
    % st=(ii-1)*sz+ii;
    % en=min(st+sz, size(rho1_mat,2));
    % Amir, try below
    % rho1_finish = rho1_mat(:,st:en);
    % rho2_finish = rho2_mat(:,st:en);
    % [hst,bn] = hist(log10(rho1_finish(:)./(rho2_finish(:))),101);

    st = ii*sz;
    rho1_finish = rho1_mat(:,st-11:st);
    rho2_finish = rho2_mat(:,st-11:st);
    
    [hst,bn] = hist(log10(rho1_finish(:)./(rho2_finish(:))),101);
    hst = hst / trapz(bn,hst);
   
    entropies(ii) = -trapz(bn, hst.*log(hst));
  
    semilogx(10.^(bn),hst,'.', 'Color', smap(ii,:));
    
    bi=st*max_batches/(sims);
 
    logr_mean_theory=bi*theory_mu/log(10);
    logr_std_theory=sqrt(bi)*theory_sigma/log(10);
    xx_theory = bn;
    % yy_theory = 1/sqrt(2*pi)/logr_std_theory*exp(-((xx_theory-logr_mean_theory).^2)/(2*(logr_std_theory)^2));
    yy_theory = 1/sqrt(2*pi)/logr_std_theory*exp(-((xx_theory-logr_mean_theory).^2)/(2*((logr_std_theory)^2)));
    % yy_theory = yy_theory / trapz(xx_theory, yy_theory);
    % semilogx(10.^[logr_mean_theory,logr_mean_theory], [0,0.3], '-', 'Color', smap(ii,:));
    semilogx(10.^(xx_theory),yy_theory,'-', 'Color', smap(ii,:), 'LineWidth', 1.5);

%     semilogx(exp(bn),-log(hst),'.-', 'Color', smap(ii,:));
    % plot_fit_gauss2(bn,hst)
end
xlim(10.^[-5,25]);
xticks(10.^[0,10,20])
set(gca,'XScale','log')
% disp(['Eff species: ' mat2str(exp(entropies))]);
% ylim([1,1e15]);
% yticks(10.^(0:5:15));
xlabel('$\rho_N/\rho_D$', 'Interpreter','latex');
ylabel('Probability', 'Interpreter','latex');
set(gca,'XScale','log');
set(gca,'YScale','linear');

c=colorbar();
c.Ticks = [0,0.5,1];
c.TickLabels={'Start', '$2\times 10^5$', '$4\times 10^5$'};
c.TickLabelInterpreter='latex';
c.Label.String='Batches';
c.Label.Color=[1,1,1];
c.Label.Position=[0.08957112187925 0.190722126321694 0];
c.Label.Units='data';
% c.Label.HorizontalAlignment='center';


print(gcf,'-dpng', '../figures/fig_stoch_histogram_noinset.png', '-r600');
print(gcf,'-dsvg', '../figures/fig_stoch_histogram_noinset.svg');




%% Drift in log space:
% figure;
s1=load('automated_runs/results/hpc_high_Omega_coexistence_10K_part1_collected_results.mat');
s1batches = s1.batches;
s1rho1_mat = s1.rho1_mat;
s1rho2_mat = s1.rho2_mat;

s2=load('automated_runs/results/hpc_high_Omega_fifth_noise_10K_precollected_results');
s2sbatches = s2.precollected_final.batch_cell{1};
s2rho1_mat = s2.precollected_final.rho1_mat;
s2rho2_mat = s2.precollected_final.rho2_mat;
load('automated_runs/results/highOmega400K_sim_group_1000_v2.mat')
str='highOmega10K-lowNoise';
det_cell=1;

col_deterministic = [17,119,51]/255;
col_weak_noise = [204,102,119]/255;
col_strong_noise = [136,34,85]/255;

%
% f=newfigure(3,2);
axes('Position',[0.52,0.55,0.21,0.32]);

% box on
set(gca,'FontSize',11)
hold on
plot((1:size(s1rho1_mat,2))*100,10.^(median(log10(s1rho1_mat./s1rho2_mat),1)),...
    '-', 'LineWidth',1.5, 'Color',col_strong_noise,...
    'DisplayName','Strong noise')
plot((1:size(s2rho1_mat,2))*100,10.^(median(log10(s2rho1_mat./s2rho2_mat),1)),...
    '-', 'LineWidth',1.5, 'Color',col_weak_noise,...
    'DisplayName','Weak noise')
semilogy((1:100:length(output_cell{det_cell}.rho(1,:))),output_cell{det_cell}.rho(1,1:100:end)./output_cell{det_cell}.rho(2,1:100:end),...
    '-', 'LineWidth',1.5, 'Color',col_deterministic,...
    'DisplayName','Deterministic')
ylim([1,1e15]);
% ylim([0,15]);
yticks(10.^[0,15]);
% xlabel('Batch', 'Interpreter','latex');
ylabel('$\rho_N/\rho_D$', 'Interpreter','latex');
set(gca,'XScale','linear');
set(gca,'YScale','log');
text(1.3e5,15,'Batches', 'Interpreter','latex')

print(gcf,'-dpng', '../figures/fig_stoch_histogram.png', '-r600');
print(gcf,'-dsvg', '../figures/fig_stoch_histogram.svg');
