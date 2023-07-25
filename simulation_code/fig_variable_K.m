%% Import and reformat data

clear;clc

%results_table = readtable(['automated_runs/results/initial_variable_K_sweep_results_table.csv'],'Delimiter',',');
results_table = readtable(['automated_runs/results/hpc_production_variable_K_sweep_results_table.csv'],'Delimiter',',');
results_table = sortrows(results_table,'id','ascend');
non_collected = strcmp(results_table.final_rho,'');
results_table.final_rho(non_collected) = {'[NaN;NaN]'};

results_table.K = cellfun(@(x) eval(x),results_table.K,'UniformOutput',false);
results_table.K2 = cell2mat(cellfun(@(x) x(2), results_table.K,'UniformOutput',false));
results_table.final_rho = cellfun(@(x) eval(x),results_table.final_rho,'UniformOutput',false);
results_table.final_rho1 = cell2mat(cellfun(@(x) x(1), results_table.final_rho,'UniformOutput',false));
results_table.final_rho2 = cell2mat(cellfun(@(x) x(2), results_table.final_rho,'UniformOutput',false));
results_table.deltaE = cellfun(@(x) eval(x),results_table.deltaE,'UniformOutput',false);
results_table.deltaE2 = cell2mat(cellfun(@(x) x(2), results_table.deltaE,'UniformOutput',false));
results_table.extinction = (results_table.final_rho1 == 0) & (results_table.final_rho2 == 0);
log10c0_range = unique(results_table.log10c0);
deltaE2_range = unique(results_table.deltaE2);
K2_range = unique(results_table.K2) - 1;

for i = 1:length(deltaE2_range)
    deltaE2 = deltaE2_range(i);
    deltaE2_ind = deltaE2 == results_table.deltaE2;
    
    for j = 1:length(log10c0_range)
        log10c0 = log10c0_range(j);
        log10c0_ind = results_table.log10c0 == log10c0;
        subtable = results_table(deltaE2_ind&log10c0_ind,:);
        
        K2_mat(i,j,:) = subtable.K2;
        S_mat(i,j,:) = subtable.end_S;
        rho1_mat(i,j,:) = subtable.final_rho1;
        rho2_mat(i,j,:) = subtable.final_rho2;
        rel_abun2_mat(i,j,:) = rho2_mat(i,j,:)./(rho1_mat(i,j,:) + rho2_mat(i,j,:));
        collected_ind = ~isnan(rel_abun2_mat(i,j,:));
        collected_S = squeeze(S_mat(i,j,collected_ind));
        collected_K2 = squeeze(K2_mat(i,j,collected_ind));
        
        if sum(collected_ind) > 0
                        
            diverse_region = find(collected_S > 0.005);
            transition_ind1 = diverse_region(1);
            transition_ind2 = diverse_region(end);
            
            transition_K2_1(i,j) = collected_K2(transition_ind1);
            transition_K2_2(i,j) = collected_K2(transition_ind2);
        end
    end
    
end

%% Get the inverse squared integrals

rho0 = 1;
for i = 1:length(deltaE2_range)
    deltaE2 = deltaE2_range(i);
    for j = 1:length(log10c0_range)
        K2 = transition_K2_2(i,j);
        log10c0 = log10c0_range(j);
        [tilde_Nr_den_sq(i,j), tilde_Nr(i,j)]  = compute_Nr_den_sq(log10c0,rho0,K2,deltaE2);
        %high_c0_Nr_den_sq_mat(i,j) = 1./((1+deltaE2).*10.^log10c0)*log(((10.^log10c0).^2)/(K2.*rho0));
    end
end
c0_range = 10.^log10c0_range;
[base_Nr_den_sq, base_Nr] = compute_Nr_den_sq(log10c0_range,rho0,1,0);
%ind = 2;
%analytical_transition_deltaK_2 = (deltaE2_range(2)./(1 + deltaE2_range(2)))...
%    .*c0_range.*log(1 + c0_range./rho0)...
%    ./log(c0_range./((transition_K2_1(ind,:)')*rho0));


%% Plot the entropy as a function of K2

colors = jet(length(log10c0_range));
linetypes = {'-','--'};
figure
hold on
for i = 1:length(deltaE2_range)
    for j = 1:length(log10c0_range)
        plot(squeeze(K2_mat(i,j,:)-1),squeeze(S_mat(i,j,:)),...
            linetypes{i},'Color',colors(j,:))
    end
end
title([linetypes{1},' is ',num2str(deltaE2_range(1)),', ',linetypes{2},' is ',num2str(deltaE2_range(2))])
set(gca,'XScale','log','YScale','log')
%set(gca,'YScale','log');
xlim([0.5e-2,1000])
ylim([1e-5,1])
xlabel('K_2-1')
ylabel('S')

colormap(jet);
h = colorbar;
h.Ticks = [0,1];
h.TickLabels = {num2str(log10c0_range(1)),num2str(log10c0_range(end))};
xlabel(h,'log10c0')

% print(gcf, '-dpng','../figures/fig_entropy_vs_K2.png','-r300');


%% Plot S phase diagram

load('automated_runs/parameters/default_params_struct.mat')

ind= 2;
newfigure(3,2);
debtor_color = [206 37 123]/255;
nondebtor_color = [15 104 194]/255;
coexistence_color = [254 97 0]/255;
hold on
c0=10.^(log10c0_range');

lower_y = 10^(-0.5);
upper_y = 1e1;

predicted_deltaK2_1 = deltaE2_range(ind)/params.E(1).*base_Nr./base_Nr_den_sq;
predicted_deltaK2_2 = deltaE2_range(ind)./(params.E(1) + deltaE2_range(ind)).*tilde_Nr(ind,:)./tilde_Nr_den_sq(ind,:);

lower_bounds = zeros(size(c0)) + lower_y;
upper_bounds = zeros(size(c0)) + upper_y;

mid_bounds1 = (transition_K2_1(ind,:)-1)./c0;
mid_bounds2 = (transition_K2_2(ind,:)-1)./c0;

x = [c0,fliplr(c0)];
y1 = [lower_bounds,fliplr(mid_bounds1)];
y2 = [mid_bounds1,fliplr(mid_bounds2)];
y3 = [mid_bounds2,fliplr(upper_bounds)];
fill(x,y2,coexistence_color,'LineWidth',1.5,'EdgeColor','k')
fill(x,y3,nondebtor_color,'LineWidth',1.5,'EdgeColor','k')
fill(x,y1,debtor_color,'LineWidth',1.5,'EdgeColor','k')

plot([1e-1,1e-1],[lower_y,upper_y],'k-','LineWidth', 1.5)
plot(c0,mid_bounds1,'k-','LineWidth', 1.5, 'DisplayName','Numerical bound 1')
plot(c0,mid_bounds2,'k-','LineWidth', 1.5,'DisplayName','Numerical bound 2')
plot(c0,predicted_deltaK2_1./c0,'w--','LineWidth', 1,'DisplayName','Analytical bound 1')
plot(c0,predicted_deltaK2_2./c0,'w--','LineWidth', 1,'DisplayName','Analytical bound 2')
% chemostat_deltaK=(deltaE2_range(ind)/params.E(1)*params.K);%./ (1-c0/(params.rho0*params.E(1)*params.tf));
% plot(c0,chemostat_deltaK./c0,'k--','LineWidth', 1,'DisplayName','Chemostat')

chemostat_deltaK=(deltaE2_range(ind)/params.E(1)*params.K)./(1-1/params.tf*log(1+c0/params.rho0));%./ (1-c0/(params.rho0*params.E(1)*params.tf));
plot(c0,chemostat_deltaK./c0,'k--','LineWidth', 1,'DisplayName','Chemostat')


% plot(c0, ones(size(c0))*deltaE2_range(ind)/params.E(1)/2,'-w');
% plot(c0, ones(size(c0))*deltaE2_range(ind)/(deltaE2_range(ind)+params.E(1))./(1+log(2)./log(c0)),'-w');


ylim([0,1500])
xlabel('$c_0/K$', 'Interpreter', 'Latex');
ylabel('$\Delta K/c_0$', 'Interpreter','latex');
set(gca,'XScale','log','YScale','log')

xlim(10.^[-1,3])
xticks(10.^[0,3])
ylim([lower_y,upper_y]);
yticks(10.^[-1,0,1]);

set(gca,'TickDir','out')
set(gca,'YminorTick','off')
set(gca,'XminorTick','off')

text(20,10.^(0.5),'Non-debtor','Interpreter','latex',...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'Color','w')
% t=text(10.^(-0.25),10.^(-0.25),'Debtor','Interpreter','latex',...
text(20,10.^(-0.4),'Debtor','Interpreter','latex',...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'Color','w')
t=text(0.65,1,'Chemostat','Interpreter','latex',...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'Color','w', 'Rotation', -58)
text(200,0.61,'Both','Interpreter','latex',...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'Color','k')

print(gcf, '-dpng','../figures/fig_variable_K.png','-r600');
print(gcf, '-dsvg','../figures/fig_variable_K.svg');

%% Plot S phase diagram

load('automated_runs/parameters/default_params_struct.mat')

ind= 2;
newfigure(3,2);
% figure(777)
debtor_color = [206 37 123]/255;
nondebtor_color = [15 104 194]/255;
coexistence_color = [254 97 0]/255;
hold on
c0=10.^(log10c0_range');

lower_y = 10^(-0.4);
upper_y = 1e1;

% predicted_deltaK2_1 = deltaE2_range(ind).*log((c0 + params.rho0)./params.rho0)./(base_Nr_den_sq./params.E(1));

predicted_deltaK2_1 = deltaE2_range(ind)/params.E(1).*base_Nr./base_Nr_den_sq;
predicted_deltaK2_2 = deltaE2_range(ind)./(params.E(1) + deltaE2_range(ind)).*tilde_Nr(ind,:)./tilde_Nr_den_sq(ind,:);

lower_bounds = zeros(size(c0)) + lower_y;
upper_bounds = zeros(size(c0)) + upper_y;

mid_bounds1 = (transition_K2_1(ind,:)-1)./c0;
mid_bounds2 = (transition_K2_2(ind,:)-1)./c0;

x = [c0,fliplr(c0)];
y1 = [lower_bounds,fliplr(mid_bounds1)];
y2 = [mid_bounds1,fliplr(mid_bounds2)];
y3 = [mid_bounds2,fliplr(upper_bounds)];
% fill(x,y2,coexistence_color,'LineWidth',1.5,'EdgeColor','k')
% fill(x,y3,nondebtor_color,'LineWidth',1.5,'EdgeColor','k')
% fill(x,y1,debtor_color,'LineWidth',1.5,'EdgeColor','k')
% 
% plot([1e-1,1e-1],[lower_y,upper_y].*c0,'k-','LineWidth', 1.5)
plot(c0,mid_bounds1.*c0,'.-k','LineWidth', 1.5, 'DisplayName','Numerical bound 1')
plot(c0,mid_bounds2.*c0,'.-k','LineWidth', 1.5,'DisplayName','Numerical bound 2')
plot(c0,predicted_deltaK2_1,'r--','LineWidth', 1,'DisplayName','Analytical bound 1')
plot(c0,predicted_deltaK2_2,'c--','LineWidth', 1,'DisplayName','Analytical bound 2')

% ylim([0,1500])
xlabel('$c_0/K$', 'Interpreter', 'Latex');
ylabel('$\Delta K$', 'Interpreter','latex');
set(gca,'XScale','log','YScale','log')

xlim(10.^[-1,3])
xticks(10.^[0,3])
% ylim([lower_y,upper_y]);
% yticks(10.^[0,1]);

set(gca,'TickDir','out')
set(gca,'YminorTick','off')
set(gca,'XminorTick','off')

% text(20,10.^(0.5),'Non-debtor','Interpreter','latex',...
%     'HorizontalAlignment','center','VerticalAlignment','middle',...
%     'Color','w')
% text(10.^(-0.25),10.^(-0.25),'Debtor','Interpreter','latex',...
%     'HorizontalAlignment','center','VerticalAlignment','middle',...
%     'Color','w')
% text(200,0.61,'Both','Interpreter','latex',...
%     'HorizontalAlignment','center','VerticalAlignment','middle',...
%     'Color','w')

% print(gcf, '-dpng','figures/fig_variable_K.png','-r600');

%% Plot Delta_K/Delta_E E/K vs. c0/K

load('automated_runs/parameters/default_params_struct.mat')

ind= 2;
newfigure(3,2);
% figure(777)
debtor_color = [206 37 123]/255;
nondebtor_color = [15 104 194]/255;
coexistence_color = [254 97 0]/255;
hold on
c0=10.^(log10c0_range');

lower_y = 10^(-0.4);
upper_y = 1e1;

% predicted_deltaK2_1 = deltaE2_range(ind).*log((c0 + params.rho0)./params.rho0)./(base_Nr_den_sq./params.E(1));

predicted_deltaK2_1 = deltaE2_range(ind)/params.E(1).*base_Nr./base_Nr_den_sq;
predicted_deltaK2_2 = deltaE2_range(ind)./(params.E(1) + deltaE2_range(ind)).*tilde_Nr(ind,:)./tilde_Nr_den_sq(ind,:);

predicted_val1 = predicted_deltaK2_1/params.K / (deltaE2_range(ind)/params.E(1));
predicted_val2 = predicted_deltaK2_2/params.K / (deltaE2_range(ind)/params.E(1));


% lower_bounds = zeros(size(c0)) + lower_y;
% upper_bounds = zeros(size(c0)) + upper_y;

mid_bounds1 = (transition_K2_1(ind,:)-params.K)/params.K / (deltaE2_range(ind)/params.E(1));
mid_bounds2 = (transition_K2_2(ind,:)-params.K)/params.K / (deltaE2_range(ind)/params.E(1));

x = [c0,fliplr(c0)];
y1 = [lower_bounds,fliplr(mid_bounds1)];
y2 = [mid_bounds1,fliplr(mid_bounds2)];
y3 = [mid_bounds2,fliplr(upper_bounds)];
% fill(x,y2,coexistence_color,'LineWidth',1.5,'EdgeColor','k')
% fill(x,y3,nondebtor_color,'LineWidth',1.5,'EdgeColor','k')
% fill(x,y1,debtor_color,'LineWidth',1.5,'EdgeColor','k')
% 
% plot([1e-1,1e-1],[lower_y,upper_y].*c0,'k-','LineWidth', 1.5)
plot(c0,mid_bounds1,'.-k','LineWidth', 1.5, 'DisplayName','Numerical bound 1')
plot(c0,mid_bounds2,'.-b','LineWidth', 1.5,'DisplayName','Numerical bound 2')
plot(c0,predicted_val1,'r--','LineWidth', 1,'DisplayName','Analytical bound 1')
plot(c0,predicted_val2,'c--','LineWidth', 1,'DisplayName','Analytical bound 2')

% ylim([0,1500])
xlabel('$c_0/K$', 'Interpreter', 'Latex');
ylabel('$\Delta_K/K / (\Delta_E/E)$', 'Interpreter','latex');
set(gca,'XScale','log','YScale','log')

% xlim(10.^[-1,3])
xticks(10.^[0,3])
% ylim([lower_y,upper_y]);
% yticks(10.^[0,1]);

set(gca,'TickDir','out')
set(gca,'YminorTick','off')
set(gca,'XminorTick','off')

% text(20,10.^(0.5),'Non-debtor','Interpreter','latex',...
%     'HorizontalAlignment','center','VerticalAlignment','middle',...
%     'Color','w')
% text(10.^(-0.25),10.^(-0.25),'Debtor','Interpreter','latex',...
%     'HorizontalAlignment','center','VerticalAlignment','middle',...
%     'Color','w')
% text(200,0.61,'Both','Interpreter','latex',...
%     'HorizontalAlignment','center','VerticalAlignment','middle',...
%     'Color','w')

% print(gcf, '-dpng','figures/fig_variable_K.png','-r600');

