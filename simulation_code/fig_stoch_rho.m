clear

  


s1=load('automated_runs/results/hpc_high_Omega_coexistence_10K_part1_collected_results.mat');
s1batches = s1.batches;
s1rho1_mat = s1.rho1_mat;
s1rho2_mat = s1.rho2_mat;

%These are the theoretical mean and std of log ratio increments ln(r_(b+1))-ln(r_b) 
%The variables used in this simulation:
c0=100 ;
rho0=1 ;
E=1 ;
deltaE=0.05 ;
Omega=0.002508737448117671 ;
max_batches=400000 ;
tf=92.0176177196183 ;
theory_mu=(Omega*tf-deltaE/E*log(1+c0/rho0));
delta_Omega2=0.000250873744811767 ;
theory_sigma=sqrt(tf^2*(delta_Omega2^2)/3);


s2=load('automated_runs/results/hpc_high_Omega_fifth_noise_10K_precollected_results');
s2sbatches = s2.precollected_final.batch_cell{1};
s2rho1_mat = s2.precollected_final.rho1_mat;
s2rho2_mat = s2.precollected_final.rho2_mat;
load('automated_runs/results/highOmega400K_sim_group_1000_v2.mat')
str='highOmega10K-lowNoise';
det_cell=1;
batch_array = (0:(length(s2rho1_mat(1,:)-2)))/length(s2rho1_mat(1,:))*max_batches ;
%These are arrays of the theoretical mean and variance in ln(r_b) over batch counts
%We used ln(r_0)=0, since I thought we had rho_N=rho_D at b=0
logr_theory_avg = batch_array*theory_mu;
logr_theory_var = batch_array*(theory_sigma^2);
xx = [-1000:1000];
frac_theory_avg = zeros(size(logr_theory_avg));
for jj=1:length(logr_theory_avg)
    mn = logr_theory_avg(jj);
    vr = logr_theory_var(jj);
    frac_theory_avg(jj) = trapz(xx, (1./(1+exp(-xx))).*exp(-((xx-mn).^2)./(2*vr))/sqrt(2*vr*pi)) ;
end

%

col_deterministic = [17,119,51]/255;
col_weak_noise = [204,102,119]/255;
col_strong_noise = [136,34,85]/255;

newfigure(3,2);
set(gca,'FontSize',11)
hold on 
plot((1:100:length(output_cell{det_cell}.rho(1,:))),output_cell{det_cell}.rho(1,1:100:end)./(output_cell{det_cell}.rho(1,1:100:end)+output_cell{det_cell}.rho(2,1:100:end)),...
    '-', 'LineWidth',1.5, 'Color',col_deterministic,...
    'DisplayName','$\Omega=\Omega^*+$ Drift')
mn = mean(s2rho1_mat./(s2rho1_mat+s2rho2_mat),1);
s = std(s2rho1_mat./(s2rho1_mat+s2rho2_mat),1)/sqrt(size(s2rho1_mat,1)-1);
x=(1:size(s2rho1_mat,2))*100;
plot(x,mn,...
    '-', 'LineWidth',1.5, 'Color', col_weak_noise,...
    'DisplayName','$\Omega\pm 10\%$ Noise')
h = plot(batch_array(10:end), frac_theory_avg(10:end), '--k', 'LineWidth',1.5, 'DisplayName','Gaussian Theory');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% h=plot(x,mn-s,...
%     '--', 'LineWidth',1.5, 'Color', col_weak_noise
%     'DisplayName','$\Omega\pm 10\%$ Noise')
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% h=plot(x,mn+s,...
%     '--', 'LineWidth',1.5, 'Color', col_weak_noise,...
%     'DisplayName','$\Omega\pm 10\%$ Noise')
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
h=patch([x fliplr(x)], [mn-s fliplr(mn+s)], col_weak_noise,'FaceAlpha',0.1,...
    'EdgeColor', col_weak_noise, 'EdgeAlpha', 0.2)
h.Annotation.LegendInformation.IconDisplayStyle = 'off';


delta_Omega2=Omega/2 ;
theory_sigma=sqrt(tf^2*(delta_Omega2^2)/3);
logr_theory_avg = batch_array*theory_mu;
logr_theory_var = batch_array*(theory_sigma^2);
xx = [-1000:1000];
frac_theory_avg = zeros(size(logr_theory_avg));
for jj=1:length(logr_theory_avg)
    mn = logr_theory_avg(jj);
    vr = logr_theory_var(jj);
    frac_theory_avg(jj) = trapz(xx, (1./(1+exp(-xx))).*exp(-((xx-mn).^2)./(2*vr))/sqrt(2*vr*pi)) ;
end

mn=mean(s1rho1_mat./(s1rho1_mat+s1rho2_mat),1);
s=std(s1rho1_mat./(s1rho1_mat+s1rho2_mat),1)/sqrt(size(s2rho1_mat,1)-1);
x=(1:size(s1rho1_mat,2))*100;
plot(x, mn,...
    '-', 'LineWidth',1.5, 'Color', col_strong_noise,...
    'DisplayName','$\Omega\pm 50\%$ Noise')
plot(batch_array, frac_theory_avg, '--k', 'LineWidth',1.5, 'DisplayName','Gaussian')

% h=plot(x, mn+s,...
%     '--', 'LineWidth',1.5, 'Color', col_strong_noise,'Alpha',0.1,...
%     'DisplayName','$\Omega\pm 50\%$ Noise')
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% h=plot(x, mn-s,...
%     '--', 'LineWidth',1.5, 'Color', col_strong_noise,'Alpha',0.1,...
%     'DisplayName','$\Omega\pm 50\%$ Noise')
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
h=patch([x fliplr(x)], [mn-s fliplr(mn+s)], col_strong_noise,'FaceAlpha',0.1,...
    'edgecolor', col_strong_noise, 'edgealpha',0.2)
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

ylim([0.5,1]);
yticks([0.5,1]);
xlabel('Batch', 'Interpreter','latex')
ylabel('Mean $\rho_N/(\rho_N+\rho_D)$', 'Interpreter','latex')
set(gca,'XScale','linear')
set(gca,'YScale','linear')
l=legend();
l.Position = [0.54 0.26 0.43 0.25];
l.FontSize=9;
l.Interpreter='latex';
print(gcf,'-dpng', '../figures/fig_stoch_rho_linear.png', '-r600');
print(gcf,'-dsvg', '../figures/fig_stoch_rho_linear.svg');


f=newfigure(3,2);
set(gca,'FontSize',11)
hold on
semilogy((1:size(s1rho1_mat,2))*100,exp(median(log(s1rho1_mat./s1rho2_mat),1)),...
    '-', 'LineWidth',1.5, 'Color',col_strong_noise,...
    'DisplayName','Strong noise')
semilogy((1:size(s2rho1_mat,2))*100,exp(median(log(s2rho1_mat./s2rho2_mat),1)),...
    '-', 'LineWidth',1.5, 'Color',col_weak_noise,...
    'DisplayName','Weak noise')
semilogy((1:100:length(output_cell{det_cell}.rho(1,:))),exp(median(log(output_cell{det_cell}.rho(1,1:100:end)./output_cell{det_cell}.rho(2,1:100:end)),1)),...
    '-', 'LineWidth',1.5, 'Color',col_deterministic,...
    'DisplayName','Deterministic')
ylim([1,1e15]);
yticks(10.^(0:5:15));
xlabel('Batch', 'Interpreter','latex');
ylabel('$\rho_N/\rho_D$', 'Interpreter','latex');
set(gca,'XScale','linear');
set(gca,'YScale','log');
l=legend();
l.Position = [0.5289 0.255 0.45 0.2847];
l.FontSize=9;
l.Interpreter='latex';
print(gcf,'-dpng', '../figures/fig_stoch_rho_log.png', '-r600');
print(gcf,'-dsvg', '../figures/fig_stoch_rho_log.svg');



