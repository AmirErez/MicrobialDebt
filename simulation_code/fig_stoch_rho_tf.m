clear

% s1=load('automated_runs/results/hpc_high_Omega_coexistence_10K_part1_collected_results.mat');
% s1batches = s1.batches;
% s1rho1_mat = s1.rho1_mat;
% s1rho2_mat = s1.rho2_mat;

s1=load('automated_runs/results/hpc_small_Omega_coexistence_tf_notlong_correctedDrift_10K_precollected_results.mat');
s1sbatches = s1.precollected_final.batch_cell{1};
s1rho1_mat = s1.precollected_final.rho1_mat;
s1rho2_mat = s1.precollected_final.rho2_mat;

s2=load('automated_runs/results/hpc_high_Omega_fifth_noise_10K_precollected_results');
s2sbatches = s2.precollected_final.batch_cell{1};
s2rho1_mat = s2.precollected_final.rho1_mat;
s2rho2_mat = s2.precollected_final.rho2_mat;
load('automated_runs/results/highOmega400K_sim_group_1000_v2.mat')
str='highOmega10K-lowNoise';
det_cell=1;

s3=load('automated_runs/results/hpc_small_Omega_coexistence_stochastic_tf_notlong_correctedDrift_10K_precollected_results.mat');
s3sbatches = s3.precollected_final.batch_cell{1};
s3rho1_mat = s3.precollected_final.rho1_mat;
s3rho2_mat = s3.precollected_final.rho2_mat;



col_deterministic = [17,119,51]/255;
col_weak_noise = [204,102,119]/255;
col_strong_noise = [136,34,85]/255;

newfigure(3,2);
set(gca,'FontSize',11)
hold on 
plot((1:size(s1rho1_mat,2))*100, mean(s1rho1_mat./(s1rho1_mat+s1rho2_mat),1),...
    '-', 'LineWidth',1.5, 'Color', col_strong_noise,...
    'DisplayName','$\Omega\pm 10\%$, $t_f\approx 9$')
plot((1:size(s2rho1_mat,2))*100,mean(s2rho1_mat./(s2rho1_mat+s2rho2_mat),1),...
    '-', 'LineWidth',1.5, 'Color', col_weak_noise,...
    'DisplayName','$\Omega\pm 10\%$, $t_f\approx 90$')
plot((1:size(s3rho1_mat,2))*100,mean(s3rho1_mat./(s3rho1_mat+s3rho2_mat),1),...
    '-', 'LineWidth',1.5, 'Color', col_weak_noise,...
    'DisplayName','$t_f\pm 10\%$, $t_f\approx 9$')
plot((1:100:length(output_cell{det_cell}.rho(1,:))),output_cell{det_cell}.rho(1,1:100:end)./(output_cell{det_cell}.rho(1,1:100:end)+output_cell{det_cell}.rho(2,1:100:end)),...
    '-', 'LineWidth',1.5, 'Color',col_deterministic,...
    'DisplayName','Deterministic')
ylim([0.5,1]);
yticks([0.5,1]);
xlim([0,3e5])
xlabel('Batch', 'Interpreter','latex')
ylabel('Mean $\rho_N/(\rho_N+\rho_D)$', 'Interpreter','latex')
set(gca,'XScale','linear')
set(gca,'YScale','linear')
l=legend();
% l.Position = [0.5289 0.255 0.4583 0.2847];
l.Location = 'southeast';
l.Interpreter='latex';
print(gcf,'-dpng', '../figures/fig_stoch_tf_rho_linear.png', '-r600');
print(gcf,'-dsvg', '../figures/fig_stoch_tf_rho_linear.svg');


f=newfigure(3,2);
set(gca,'FontSize',11)
hold on
semilogy((1:size(s1rho1_mat,2))*100,exp(median(log(s1rho1_mat./s1rho2_mat),1)),...
    '-', 'LineWidth',1.5, 'Color',col_strong_noise,...
    'DisplayName','$\Omega\pm 10\%$, $t_f\approx 9$')
semilogy((1:size(s2rho1_mat,2))*100,exp(median(log(s2rho1_mat./s2rho2_mat),1)),...
    '-', 'LineWidth',1.5, 'Color',col_weak_noise,...
    'DisplayName','$\Omega\pm 10\%$, $t_f\approx 90$')
semilogy((1:size(s3rho1_mat,2))*100,exp(median(log(s3rho1_mat./s3rho2_mat),1)),...
    '-', 'LineWidth',1.5, 'Color',col_weak_noise,...
    'DisplayName','$t_f\pm 10\%$, $t_f\approx 9$')
semilogy((1:100:length(output_cell{det_cell}.rho(1,:))),exp(median(log(output_cell{det_cell}.rho(1,1:100:end)./output_cell{det_cell}.rho(2,1:100:end)),1)),...
    '-', 'LineWidth',1.5, 'Color',col_deterministic,...
    'DisplayName','Deterministic')
ylim([1,1e10]);
xlim([0,3e5])
yticks(10.^(0:5:15));
xlabel('Batch', 'Interpreter','latex');
ylabel('$\rho_N/\rho_D$', 'Interpreter','latex');
set(gca,'XScale','linear');
set(gca,'YScale','log');
l=legend();
% l.Position = [0.5289 0.255 0.4583 0.2847];
l.Position = [0.078703703703704,0.624999999999999,0.535601843299372,0.373690748712962];
l.Interpreter='latex';
print(gcf,'-dpng', '../figures/fig_stoch_tf_rho_log.png', '-r600');
print(gcf,'-dsvg', '../figures/fig_stoch_tf_rho_log.svg');


