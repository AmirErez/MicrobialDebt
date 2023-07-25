%% Run new numerical code

clear;clc

params.p = 1;
params.m = 2;
params.K = 1;
params.max_batches = 1e6;
%params.Q =
params.E = [1;1];
params.P = 1;
params.d = [0;0];
params.rho0 = 1;
params.rho_at_t0 = [1/2;1/2];
params.psi = [1;1.05];
params.theta = [0;0];
params.log10c0 = 2;
params.alpha = [1;1];
params.errtype = 1;
params.dilution_method = 'floating';
plt = 0;

%center = 0.09105;
%del = 0.0005;
%Omega_vec = linspace(0.09105-del,0.09105+del,40);
Omega_vec = linspace(0.875,0.9,40);
%c0_vec = [-2,-1,0,1];
%log10c0_vec = linspace(-3,-4,2);
%log10c0_vec = [-3,-2,-1,0];
log10c0_vec = 1;
%log10c0_vec = [0.05,0.1,0.15,0.25];

tic
for j = 1:length(log10c0_vec)
    params.log10c0 = log10c0_vec(j);
    
    for i = 1:length(Omega_vec)
        params.Omega = [1;Omega_vec(i)];
        tic
        output = serialdil_odesolver(params, plt);
        rho_vec(:,i,j) = output.rho(:,end);
        disp(['Completed run ',num2str(i),' with log10(c_0) = ',...
            num2str(params.log10c0),' in ',num2str(toc,2)])
    end
end
toc

%save('initial_phase_diagram_deltaE_2_low_c0.mat')
%save('initial_phase_diagram_deltaE_2_high_c0.mat')
save('sim_data/new_code_wide_check_3.mat')


%%
H = @(x) nansum(-x.*log2(x));

for j = 1:length(log10c0_vec)
    for i = 1:length(Omega_vec)
        H_mat(i,j) = H(rho_vec(:,i,j)./sum(rho_vec(:,i,j)));
        
    end
end
H_mat = fliplr(H_mat);
log10c0_vec = fliplr(log10c0_vec);

figure
hold on
colors = jet(length(log10c0_vec));
for i = 1:1:length(log10c0_vec)
    %plot(c0_vec,H_mat(i,:),'-','Color',colors(i,:),...
    %    'DisplayName',num2str(dilution_2_vec(i)));
    plot(Omega_vec,H_mat(:,i),'.-','Color',colors(i,:),...
        'DisplayName',num2str(log10c0_vec(i)),'LineWidth',1);
    %'DisplayName',num2str(log2((10.^c0_vec(i) + params.rho0)./params.rho0)));
end
%xlim([1,50])
%ylim([1e-7,1])
legend()
ylabel('H')
xlabel('d_2/d_1')

%% PLOT

clear;clc
load('new_code_wide_check_2.mat')

H = @(x) nansum(-x.*log2(x));

for j = 1:length(log10c0_vec)
    for i = 1:length(Omega_vec)
        H_mat(i,j) = H(rho_vec(:,i,j)./sum(rho_vec(:,i,j)));
        
    end
end
H_mat = fliplr(H_mat);
log10c0_vec = fliplr(log10c0_vec);

H_cutoff = 1e-3;
figure
%imagesc(log((10.^c0_vec+params.rho0)./params.rho0),dilution_2_vec-1,log10(H_mat))
imagesc(log10c0_vec,Omega_vec./1,log10(H_mat))
set(gca,'YDir','normal')
%xticks([1,size(H_mat,2)])
%xticklabels({'-4','-3'})
xlabel('')
%yticks([1,size(H_mat,1)])
%yticklabels({'10','1'})
ylabel('d_2/d_1')


%% plot lines

H = @(x) nansum(-x.*log2(x));

figure
hold on
colors = jet(length(mult_2_vec));
for i = 1:1:length(mult_2_vec)
    %plot(c0_vec,H_mat(i,:),'-','Color',colors(i,:),...
    %    'DisplayName',num2str(dilution_2_vec(i)));
    plot(log2((10.^log10c0_vec + params.rho0)./params.rho0),H_mat(i,:),'-','Color',colors(i,:),...
        'DisplayName',num2str(mult_2_vec(i)));
end
legend()
ylabel('H')
xlabel('Generations')
%set(gca,'YScale','linear')

%% plot lines in dilution space
clear;clc
load('new_code_wide_check_2.mat')

H = @(x) nansum(-x.*log2(x));

for j = 1:length(log10c0_vec)
    for i = 1:length(Omega_vec)
        H_mat(i,j) = H(rho_vec(:,i,j)./sum(rho_vec(:,i,j)));
        
    end
end
H_mat = fliplr(H_mat);
log10c0_vec = fliplr(log10c0_vec);

figure
hold on
colors = jet(length(log10c0_vec));
for i = 1:1:length(log10c0_vec)
    %plot(c0_vec,H_mat(i,:),'-','Color',colors(i,:),...
    %    'DisplayName',num2str(dilution_2_vec(i)));
    plot(Omega_vec,H_mat(:,i),'.-','Color',colors(i,:),...
        'DisplayName',num2str(log10c0_vec(i)),'LineWidth',1);
    %'DisplayName',num2str(log2((10.^c0_vec(i) + params.rho0)./params.rho0)));
end
%xlim([1,50])
%ylim([1e-7,1])
legend()
ylabel('H')
xlabel('d_2/d_1')
%set(gca,'YScale','log')
%set(gca,'XScale','log')
%set(gca,'YScale','linear')

%% Check low c0 scaling at specified d

for i = 1:length(c0_vec)
    prediction_d = exp(2*(10.^log10c0_vec(i)));
    pred_H(i) = interp1(mult_2_vec,H_mat(:,i),prediction_d);
end

figure
loglog(10.^log10c0_vec,pred_H,'k-')


%% Replot low c0 limit (oh god send help)

clear;clc

H = @(x) nansum(-x.*log2(x));

%load('initial_phase_diagram_deltaE_2_low_c0.mat')
%load('initial_phase_diagram_deltaE_2_high_c0.mat')
load('initial_phase_diagram_deltaE_2_mid_c0.mat')

colors = jet(length(log10c0_vec));
figure
hold on
for j = 1:length(log10c0_vec)
    %     params.log10c0 = log10c0_vec(j);
    predict_d(j) = exp(-params.delta_E(2)*10.^(log10c0_vec(j)));
    %     mult_2_vec{j} = mult_2*10.^(linspace(-1,1,21));
    for i = 1:length(mult_2_vec{j})
        
        H_mat(i,j) = H(rho_vec(:,i,j));
        
        %         params.dilution_factor = [1;mult_2_vec{j}(i)];
        %         output = serialdil_odesolver(params, plt);
        %         rho_vec(:,i,j) = output.rho(:,end);
        %         disp(['Completed run ',num2str(i),' with log10(c_0) = ', num2str(params.log10c0)])
    end
    %     plot(mult_2_vec{j}/predict_d(j),H_mat(:,j),'.-','Color',colors(j,:),'DisplayName',num2str(log10c0_vec(j)))
    %     predict_H(j) = interp1(mult_2_vec{j},H_mat(:,j), predict_d(j));
    %     plot(predict_d(j)/predict_d(j), predict_H(j),'*','MarkerSize', 15, 'Color', colors(j,:))
    
    plot(mult_2_vec{j},H_mat(:,j),'.-','Color',colors(j,:),'DisplayName',num2str(log10c0_vec(j)))
    predict_H(j) = interp1(mult_2_vec{j},H_mat(:,j), predict_d(j));
    plot(predict_d(j), predict_H(j),'*','MarkerSize', 15, 'Color', colors(j,:))
    
    
end
legend()
set(gca,'XScale','log')
set(gca,'YScale','log')


%% Look at the stochastic code
clear;clc
params.p = 1;
params.m = 2;
params.K = 1;
params.tf = compute_single_species_tf(0.99,4,1,1,'hill',1);
params.E = [1;1];
params.P = 1;
params.rho0 = 1;
params.rho_at_t0 = [0.5*params.rho0;0.5*params.rho0];
params.deltaE = [0;1];
params.log10c0 = -2.368421052631580;
params.Omega = [0;0];
params.omega = [0;0];
params.alpha = [1;1];
params.errtype = 1;
params.dilution_method = 'floating';
params.gtype = 'hill';
params.h = 1;
%params.environment_type = 'deterministic';
%params.environment_type = 'stochastic-c0';
params.environment_type = 'stochastic-tf';
%params.log10c0_var = 2;
params.tf_var = 1;
plt = 1;

omega2_range = linspace(0.0001,0.001,40);
enviro_type_cell = {'deterministic','stochastic-tf'};
for i = 9:9%length(omega2_range)
    for j = 1:length(enviro_type_cell)
        params.omega = [0;omega2_range(i)];
        params.environment_type = enviro_type_cell{j};
        if j == 1
            params.max_batches = 1e6;
        elseif j == 2
            params.max_batches = round(30*length(output_cell{i,1}.t));
        end
        output_cell{i,j} = serialdil_odesolver(params, plt);
        rho2_rel_abun(i,j) = output_cell{i,j}.rho(2,end)/sum(output_cell{i,j}.rho(:,end));
        disp([enviro_type_cell{j},' ',num2str(i),' completed.'])
    end
end


%%
clear;clc
n = 100;
for i = 1:n
    [output_cell{i},params_cell{i}] = run_random_2nutrient_community(50,1,1,'punctuated','deterministic');
    end_S(i) = output_cell{i}.ShannonS(end);
end

max(end_S)


%%
load('automated_runs/parameters/default_params_struct.mat')

params.m = 2;
params.log10c0 = 1;
params.rho0 = 1;
params.rho_at_t0 = zeros(params.m,1) + params.rho0/params.m;
params.E = ones(params.m,1);

%params.deltaE = [0.912217883562916;0.479489268503963];
%params.K = [1.870600909001257; 0.328622767363803];
params.deltaE = [0.91;0.48];
params.K = [1.9; 0.35];
params.p = 1;
params.alpha = ones(params.m,1);
params.P = 1;

params.environment_type = 'deterministic';
params.max_batches = 100000;

plt = 1;
output = serialdil_odesolver(params,plt);

%% Look at integral dependence on tfinal

clear;clc
params.p = 1;
params.m = 1;
params.K = 1;
params.mean_tf = compute_single_species_tf(0.99,4,1,1,'hill',1);
params.E = 1;
params.P = 1;
params.rho0 = 1;
params.rho_at_t0 = params.rho0;
params.deltaE = 0;
params.Omega = 0;
params.omega = 0;
params.alpha = 1;
params.errtype = 1;
params.dilution_method = 'floating';
params.gtype = 'hill';
params.h = 1;
params.environment_type = 'deterministic';
log10c0_vec = linspace(-1,3,10);
tf_vec = logspace(log10(0.001*params.mean_tf),log10(1000*params.mean_tf),0.25e3);
plt = 0;

params.tf_var = 1;
logn_mode_tf = log(params.mean_tf) - (params.tf_var^2)/2;
tf_pdf = lognpdf(tf_vec,logn_mode_tf,params.tf_var);

for j = 1:length(log10c0_vec)
    params.log10c0 = log10c0_vec(j);
    for i = 1:length(tf_vec)
        params.tf = tf_vec(i);
        [rho_sigma{i}, Nr(i), c_i{i}, t{i}] = multispeciesbatch_odesolver(params,plt);
        added_rho(i) = 10.^params.log10c0 - c_i{i}(end);
        final_rho(i) = params.rho0 + added_rho(i);
    end
    
    growth_minus_dilution(j) = trapz(tf_vec,tf_pdf.*added_rho) - (10.^params.log10c0);
end

figure
hold on
%plot(tf_vec, exp(-params.omega.*tf_vec),'k--')
plot(tf_vec,exp_Nr,'r--')
set(gca,'XScale','log')
%plot(tf_vec,exp(Nr - params.omega.*tf_vec),'r-')

%% Check the low deltaK relation

clear;clc
params.p = 1;
params.m = 1;
params.K = 1;
params.tf = compute_single_species_tf(0.99,4,1,1,'hill',1);
params.E = 1;
params.P = 1;
params.rho0 = 1;
params.rho_at_t0 = params.rho0;
params.deltaE = 0;
params.log10c0 = 1;
params.Omega = 0;
params.omega = 0;
params.alpha = 1;
params.errtype = 1;
params.dilution_method = 'floating';
params.gtype = 'hill';
params.h = 1;
params.environment_type = 'deterministic';
plt = 0;
[rho_sigma, Nr, c_i, t] = multispeciesbatch_odesolver(params,plt);

gfun = @(x,K) x./(x + K).^2;
Nr_sq = trapz(t,transpose(gfun(c_i,params.K))); %Compute integral of growth

deltaK = 0.1;
deltaE = params.E*deltaK*Nr_sq/log((10.^params.log10c0 + params.rho0)/params.rho0);


% now check for coexistence
params.p = 1;
params.m = 2;
params.max_batches = 1e6;
params.K = [1; 1+ 0.99*deltaK];
params.tf = compute_single_species_tf(0.99,4,1,1,'hill',1);
params.E = [1;1];
params.P = 1;
params.rho0 = 1;
params.rho_at_t0 = 0.5*[params.rho0;params.rho0];
params.deltaE = [0;deltaE];
params.log10c0 = 1;
params.Omega = [0;0];
params.omega = [0;0];
params.alpha = [1;1];
params.errtype = 1;
params.dilution_method = 'floating';
params.gtype = 'hill';
params.h = 1;
params.environment_type = 'deterministic';

plt = 1;
output = serialdil_odesolver(params,plt);

%% Look for quasi-steady-state in run 6502


group_num = 6502;
table_file = 'automated_runs/parameters/initial_continuous_stochastic_sweep.csv';
run_per_job = 1;
%Load parameter table
opts = detectImportOptions(table_file);
opts.VariableTypes(1:end-1) = {'char'};  
%opts.
params_table = readtable(table_file,opts,'ReadRowNames',false,...
    'ReadVariableNames',true);

%Initialize storage variables
output_cell = {};

%Compute start and end simulation indices
start_end(1) = (group_num-1)*run_per_job + 1; 
start_end(2) = start_end(1) + run_per_job - 1;

disp(['Running simulations ',num2str(start_end(1)), ' to ', num2str(start_end(2))]);

plt = 0;

i = 6502;
%Read in default parameter structure
load('parameters/default_params_struct.mat');

%Replace parameters
non_id_vars = params_table.Properties.VariableNames(1:(end-1));
matching_row = params_table.id == i;
for j = 1:length(non_id_vars)
    matching_entry = params_table{matching_row,non_id_vars{j}};
    if contains(non_id_vars{j},{'dilution_method','environment_type','gtype'})
        params.(non_id_vars{j}) = matching_entry{1};
    else
        params.(non_id_vars{j}) = eval(matching_entry{1});
    end
end

params.max_batches = 3*params.max_batches;

disp(['Running simulation ',num2str(i),'.'])

%Run simulation
for i = 1:10
    output_cell{i} = serialdil_odesolver(params,0);
end


%% Plot the 

colors = jet(10);
for i = 1:10
    figure
hold on
    plot(output_cell{i}.rho(1,:),'-','Color',colors(1,:),'DisplayName','non-debtor')
    plot(output_cell{i}.rho(2,:),'--','Color',colors(10,:),'DisplayName','debtor')
set(gca,'XScale','log','YScale','log')
legend()
end


%% Check equivalency condition for high c0 at boundary
clear;clc
load('automated_runs/parameters/default_params_struct.mat')

params.deltaE = [0;0.05];
params.log10c0 = 2;
%params.tf = 10*params.tf;
params.tf = 9.2018; 
params.log10c0_var = 1;
params.rho_at_t0 = [1-0.5;0.5];
%params.Omega = [0;params.deltaE(2)/params.tf*log(1+(10.^params.log10c0)/params.rho0)];
params.Omega = [0;1e-6 + (params.deltaE(2)./params.E(1)).*log(1+(10.^params.log10c0))/(params.rho0.*params.tf)];
%params.delta_Omega2 = 0.5*params.Omega(2);
params.delta_tf = 0.5*params.tf; 
%params.environment_type = 'deterministic';
%params.environment_type = 'stochastic-c0';
%params.environment_type = 'stochastic-Omega';
params.max_batches = 100000;
%states = {'deterministic','stochastic-Omega','stochastic-Omega-bernoulli'};
states = {'deterministic','stochastic-tf'};
for i = 1:2
    params.environment_type = states{i};
    output_cell{i} = serialdil_odesolver(params, 1);
end
