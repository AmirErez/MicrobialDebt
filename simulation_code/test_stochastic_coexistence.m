clear all
% load('automated_runs/results/hpc_check_stochastic_coexistence_collected_results.mat')
% load('automated_runs/results/highOmega_sim_group_1000.mat')
% str='high-omega';
% det_cell=1;

% load('automated_runs/results/hpc_bernoulli_stochastic_coexistence_collected_results.mat')
% load('automated_runs/results/highOmega_sim_group_1000.mat')
% str='bernoulli';
% det_cell=1;

% load('automated_runs/results/hpc_critical_stochastic_coexistence_collected_results.mat')
% load('automated_runs/results/critical_sim_group_1000.mat')
% str='critical';
% det_cell=1;

% load('automated_runs/results/hpc_opposite_stochastic_coexistence_collected_results.mat')
% load('automated_runs/results/opposite_sim_group_1000.mat')
% str='opposite';
% output_cell{det_cell}.rho1 = output_cell{det_cell}.rho(1,:)
% output_cell{det_cell}.rho2 = output_cell{det_cell}.rho(2,:)
% det_cell=2;

% s1=load('automated_runs/results/hpc_high_Omega_coexistence_10K_part1_collected_results.mat');
% batches = s1.batches;
% % rho1_mat = s1.rho1_mat;
% % rho2_mat = s1.rho2_mat;
% s2=load('automated_runs/results/hpc_high_Omega_coexistence_10K_part2_precollected_results.mat');
% % rho1_mat = s2.precollected_final.rho1_mat;
% % rho2_mat = s2.precollected_final.rho2_mat;
% rho1_mat = [s1.rho1_mat;s2.precollected_final.rho1_mat];
% rho2_mat = [s1.rho2_mat;s2.precollected_final.rho2_mat];
% load('automated_runs/results/highOmega_sim_group_1000.mat')
% str='highOmega10K';
% det_cell=1;


% s1=load('automated_runs/results/hpc_high_Omega_fifth_noise_10K_precollected_results');
% batches = s1.precollected_final.batch_cell{1};
% rho1_mat = s1.precollected_final.rho1_mat;
% rho2_mat = s1.precollected_final.rho2_mat;
% load('automated_runs/results/highOmega_sim_group_1000.mat')
% str='highOmega10K-lowNoise';
% det_cell=1;

% load('automated_runs/results/hpc_moderate_Omega_stochastic_collected_results.mat');
% load('automated_runs/results/highOmega_sim_group_1000.mat')
% str='moderateOmega_small';
% det_cell=1;


s1=load('automated_runs/results/hpc_tf_noise_deviation_check_precollected_results.mat');
batches = s1.precollected_final.batch_cell{1};
rho1_mat = s1.precollected_final.rho1_mat;
rho2_mat = s1.precollected_final.rho2_mat;
load('automated_runs/results/sim_group_hpc_tf_noise_deviation_check.mat')
str='tf_noise_as_highOmega';
det_cell=1;


% Debtor starts small, noisy t_f
% s1=load('automated_runs/results/hpc_tf_noise_different_starting_precollected_results.mat');
% batches = s1.precollected_final.batch_cell{1};
% rho1_mat = s1.precollected_final.rho1_mat;
% rho2_mat = s1.precollected_final.rho2_mat;
% load('automated_runs/results/hpc_tf_noise_different_starting-deterministic.mat')
% str='tf_noise_as_highOmega-different-starting';
% det_cell=1;

% % Different initial abundances
% % load('automated_runs/results/hpc_deterministic_stochastic_IC_sweep_collected_results.mat')
% % rho1_mat = rho1_mat(1:1000,:); % Stochastic
% % rho2_mat = rho2_mat(1:1000,:);
% % rho1_mat = rho1_mat(1001:end,:);  % Deterministic
% % rho2_mat = rho2_mat(1001:end,:);
% % load('automated_runs/results/highOmega_sim_group_1000.mat')
% % str='Initial-abundances';
% % det_cell=1;

f=figure(); hold on 
% errorbar(mean((rho1_mat./(rho1_mat+rho2_mat)),1),std((rho1_mat./(rho1_mat+rho2_mat)),1),'.-')
f.Name=str;
subplot(1,2,1)
plot(mean(rho1_mat./(rho1_mat+rho2_mat),1),'.-')
ylabel('$E[\rho_1/(\rho_1+\rho_2)]$','Interpreter','latex')
hold on
% plot(output_cell{det_cell}.rho(1,1:100:end)./(output_cell{det_cell}.rho(1,1:100:end)+output_cell{det_cell}.rho(2,1:100:end)),'.-')
plot(precollected_cell{det_cell}.rho1./(precollected_cell{det_cell}.rho1+precollected_cell{det_cell}.rho2),'.-')

subplot(1,2,2);
semilogy(exp(median(log(rho1_mat./rho2_mat),1)),'.-')
hold on
ylabel('Median $\ln(\rho_1/\rho_2)$', 'Interpreter','latex')
% plot(exp(median(log(  output_cell{det_cell}.rho1 ./ output_cell{det_cell}.rho2 ),1)),'.-')
plot(exp(median(log(  precollected_cell{det_cell}.rho1 ./ precollected_cell{det_cell}.rho2 ),1)),'.-')

% plot(median(log(output_cell{det_cell}.rho(1,1:100:end)./output_cell{det_cell}.rho(2,1:100:end)),1),'.-')
% ylim([1e-15,1e15]);
% ylim([0,10])
set(gca,'XScale','linear')
set(gca,'YScale','log')
hold on

%
% sz=400;
% nsteps=10;
sz=500;
nsteps=8;
smap=colormap(copper(nsteps));
sz=size(rho1_mat,2)/nsteps;
figure(222);
for ii=1:nsteps
    st=(ii-1)*sz+ii;
    en=min(st+sz, size(rho1_mat,2));
    rho1_finish = rho1_mat(:,st:en);
    rho2_finish = rho2_mat(:,st:en);
    
    [hst,bn] = hist(log(rho1_finish(:)./(rho2_finish(:))),101);
    hst = hst / trapz(bn,hst);
    if(ii==nsteps)
%         semilogx(exp(bn),hst,'.-', 'Color',smap(ii,:), 'LineWidth',2, 'DisplayName',str);
    end
    hold on
    disp(['Entropy: ' num2str(-trapz(bn, hst.*log(hst)))]);
    semilogx(exp(bn),hst,'.-', 'Color', smap(ii,:));
%     semilogx(exp(bn),-log(hst),'.-', 'Color', smap(ii,:));
    % plot_fit_gauss2(bn,hst)
end
xlim(10.^[-40,40]);
set(gca,'XScale','log')
% legend()
% set(gca,'YScale','log')

%% Mean square displacement

% sz=size(rho1_mat,2);

rho1_finish = rho1_mat(:,st:en);
rho2_finish = rho2_mat(:,st:en);

logratio = log(rho1_mat./rho2_mat);
mean_logratio = mean(logratio,1);
x_minus_mean = logratio-logratio(:,1)-mean_logratio;
sq_x_minus_mean = x_minus_mean.^2;

mean_sq_x_minus_mean = mean(sq_x_minus_mean,1);
f = figure(111);
f.Name = 'Mean (x(t)-x(0)-E[x(t)])^2 with x= log \rho_1/\rho_2';
subplot(3,1,1);
plot(sqrt(mean_sq_x_minus_mean),'.-')
ylabel('$\langle (x(t)-x_0-E[x(t)])^2 \rangle^{1/2}$', 'Interpreter','latex')
subplot(3,1,2);     
plot(mean_logratio,'.-')
ylabel('$\langle \log{\frac{\rho_1}{\rho_2}} \rangle$', 'Interpreter','latex')

subplot(3,1,3);     
plot(mean_sq_x_minus_mean./mean_logratio,'.-')
ylabel('Var/mean', 'Interpreter','latex')

%% Plot individual convergences of the deterministic at different initial abundances
load('automated_runs/results/hpc_deterministic_stochastic_IC_sweep_collected_results.mat')
rho1_mat = rho1_mat(1001:end,:); % Deterministic only
rho2_mat = rho2_mat(1001:end,:);
f=figure(); hold on 
% errorbar(mean((rho1_mat./(rho1_mat+rho2_mat)),1),std((rho1_mat./(rho1_mat+rho2_mat)),1),'.-')
f.Name='Individual Determinstic varying initial conditions';
cmap = colormap(parula(size(rho1_mat,1)));
for ii=1:size(rho1_mat,1)
    semilogy(rho1_mat(ii,:)./rho2_mat(ii,:),'.-', 'Color', cmap(ii,:))
    hold on
end
ylim([1e-30,1e30]);
set(gca,'XScale','linear')
set(gca,'YScale','log')
hold on

