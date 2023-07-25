
%% Set up for single parameter sweep with continuous death

clear;clc

n_ind1 = 20;
n_ind2 = 2;
n_ind3 = 300;
n_ind4 = 2;

table_columns = {'log10c0','Omega','deltaE','gtype','h'};
n_par = length(table_columns);
n_tot = n_ind1*n_ind2*n_ind3*n_ind4;

load('automated_runs/parameters/default_params_struct.mat');

parameter_table = cell(n_tot,n_par);

parameter_table = cell2table(parameter_table,'VariableNames',table_columns);

log10c0_range = linspace(-3,3,n_ind1);
deltaE2_range = [0.05,1];
gtype_range = {'linear','hill'};
h_range = [NaN,1];

tn = 1;
for i = 1:n_ind1
    log10c0 = log10c0_range(i);
    for k = 1:n_ind2
        deltaE2 = deltaE2_range(k);
        high_c0_Omega = deltaE2*log(10.^log10c0)/params.tf;
        low_c0_Omega = deltaE2*(10.^log10c0)/params.tf;
        
        if log10c0 <= 0
            predicted_Omega = low_c0_Omega;
        else
            predicted_Omega = high_c0_Omega;
        end
        
        min_Omega = 0.1*predicted_Omega;
        max_Omega = min([10*predicted_Omega,1+deltaE2]);
        Omega_2_range = logspace(log10(min_Omega),log10(max_Omega),n_ind3);
        
        for j = 1:n_ind3
            for w = 1:n_ind4
                parameter_table{tn,'log10c0'} = {num2str(log10c0,16)};
                parameter_table{tn,'Omega'} ={ mat2str([0;Omega_2_range(j)],16)};
                parameter_table{tn,'deltaE'} = { mat2str([0;deltaE2],16)};
                parameter_table{tn,'gtype'} = gtype_range(w);
                parameter_table{tn,'h'} = {num2str(h_range(w))};
                tn = tn + 1;
            end
        end
    end
end

parameter_table{:,'id'} = (1:n_tot)';

writetable(parameter_table,'automated_runs/parameters/initial_punctuated_sweep.csv','Delimiter',',')

%% Set up for single parameter sweep with continuous death

clear;clc

n_ind1 = 20;
n_ind2 = 2;
n_ind3 = 300;
n_ind4 = 2;

table_columns = {'log10c0','omega','deltaE','gtype','h'};
n_par = length(table_columns);
n_tot = n_ind1*n_ind2*n_ind3*n_ind4;

load('automated_runs/parameters/default_params_struct.mat');

parameter_table = cell(n_tot,n_par);

parameter_table = cell2table(parameter_table,'VariableNames',table_columns);

log10c0_range = linspace(-3,3,n_ind1);
deltaE2_range = [0.05,1];
gtype_range = {'linear','hill'};
h_range = [NaN,1];

tn = 1;
for i = 1:n_ind1
    log10c0 = log10c0_range(i);
    for k = 1:n_ind2
        deltaE2 = deltaE2_range(k);
        high_c0_omega = deltaE2*log(10.^log10c0)/params.tf;
        low_c0_omega = deltaE2*(10.^log10c0)/params.tf;
        
        if log10c0 <= 0
            predicted_omega = low_c0_omega;
        else
            predicted_omega = high_c0_omega;
        end
        
        min_omega = 0.1*predicted_omega;
        max_omega = min([10*predicted_omega,1+deltaE2]);
        omega_2_range = logspace(log10(min_omega),log10(max_omega),n_ind3);
        
        for j = 1:n_ind3
            for w = 1:n_ind4
                parameter_table{tn,'log10c0'} = {num2str(log10c0,16)};
                parameter_table{tn,'omega'} ={ mat2str([0;omega_2_range(j)],16)};
                parameter_table{tn,'deltaE'} = { mat2str([0;deltaE2],16)};
                parameter_table{tn,'gtype'} = gtype_range(w);
                parameter_table{tn,'h'} = {num2str(h_range(w))};
                tn = tn + 1;
            end
        end
    end
end


parameter_table{:,'id'} = (1:n_tot)';

writetable(parameter_table,'automated_runs/parameters/initial_continuous_sweep.csv','Delimiter',',')

%% Stochastic punctuated dilution sweep with stochastic c0

clear;clc

parameter_table = readtable(['automated_runs/results/initial_punctuated_sweep_results_table.csv']);
parameter_table = parameter_table(strcmp(parameter_table.gtype,'hill'),:);
parameter_table.Properties.VariableNames{8} = 'max_batches';
parameter_table.max_batches(isnan(parameter_table.max_batches)) = 1e5;
parameter_table.max_batches = round(parameter_table.max_batches*6);
parameter_table.tf_var = zeros(size(parameter_table,1),1) + 1;
environment_type = cell(size(parameter_table,1),1);
environment_type(:) = {'stochastic-c0'};
parameter_table.environment_type = environment_type;
parameter_table = parameter_table(:,[1,2,3,4,5,8,11,12,6]);
parameter_table{:,'id'} = (1:size(parameter_table,1))';

writetable(parameter_table,'automated_runs/parameters/initial_punctuated_stochastic_sweep.csv','Delimiter',',')


%% Stochastic continuous dilution sweep with stochastic tf

clear;clc

parameter_table = readtable(['automated_runs/results/initial_continuous_sweep_results_table.csv']);
parameter_table = parameter_table(strcmp(parameter_table.gtype,'hill'),:);
parameter_table.Properties.VariableNames{8} = 'max_batches';
parameter_table.max_batches(isnan(parameter_table.max_batches)) = 1e5;
parameter_table.max_batches = round(parameter_table.max_batches*10);
parameter_table.max_batches(parameter_table.max_batches < 2000) = 2000;
parameter_table.tf_var = zeros(size(parameter_table,1),1) + 1;
environment_type = cell(size(parameter_table,1),1);
environment_type(:) = {'stochastic-tf'};
parameter_table.environment_type = environment_type;
parameter_table = parameter_table(:,[1,2,3,4,5,8,11,12,6]);
parameter_table{:,'id'} = (1:size(parameter_table,1))';

writetable(parameter_table,'automated_runs/parameters/initial_continuous_stochastic_sweep.csv','Delimiter',',')


%% Set up for initial variable K sweep

clear;clc

n_ind1 = 40;
n_ind2 = 2;
n_ind3 = 300;

table_columns = {'log10c0','K','deltaE'};
n_par = length(table_columns);
n_tot = n_ind1*n_ind2*n_ind3;

load('automated_runs/parameters/default_params_struct.mat');

parameter_table = cell(n_tot,n_par);

parameter_table = cell2table(parameter_table,'VariableNames',table_columns);

log10c0_range = linspace(-3,3,n_ind1);
deltaE2_range = [0.05,1];

Nr_den_sq = compute_Nr_den_sq(log10c0_range,params.rho0);

tn = 1;
for i = 1:n_ind1
    log10c0 = log10c0_range(i);
    for k = 1:n_ind2
        deltaE2 = deltaE2_range(k);
        
        predicted_deltaK = ...
            deltaE2_range(k).*log((10.^log10c0_range(i) + params.rho0)./params.rho0)...
            ./(Nr_den_sq(i)./params.E(1));

        K2_range = linspace(1+0.1*predicted_deltaK,1+6*predicted_deltaK,n_ind3);
        
        for j = 1:n_ind3
            parameter_table{tn,'log10c0'} = {num2str(log10c0,16)};
            parameter_table{tn,'K'} ={ mat2str([1;K2_range(j)],16)};
            parameter_table{tn,'deltaE'} = { mat2str([0;deltaE2],16)};
            
            tn = tn + 1;
        end
    end
end

parameter_table{:,'id'} = (1:n_tot)';

writetable(parameter_table,'automated_runs/parameters/initial_variable_K_sweep.csv','Delimiter',',')

%% Checking for true stochastic coexistence

clear;clc
load('automated_runs/parameters/default_params_struct.mat');
table_columns = {'log10c0','deltaE','Omega','delta_Omega2','max_batches','environment_type','tf',};
n_par = length(table_columns);
n_tot = 1000;
params.deltaE = [0;0.05];
params.log10c0 = 2;
params.tf = 10*params.tf;
params.Omega = [0;1e-6 + (params.deltaE(2)./params.E(1)).*log(1+(10.^params.log10c0))/(params.rho0.*params.tf)];
params.delta_Omega2 = 0.5*params.Omega(2);
params.environment_type = 'stochastic-Omega';
params.max_batches = 2e5;

parameter_table = cell(n_tot,n_par);

parameter_table = cell2table(parameter_table,'VariableNames',table_columns);

for i = 1:(n_tot-1)
    
    parameter_table{i,'log10c0'} = {num2str(params.log10c0,16)};
    parameter_table{i,'Omega'} ={ mat2str(params.Omega,16)};
    parameter_table{i,'deltaE'} = { mat2str(params.deltaE,16)};
    parameter_table{i,'delta_Omega2'} = { mat2str(params.delta_Omega2)};
    parameter_table{i,'max_batches'} = { mat2str(params.max_batches)};
    parameter_table{i,'environment_type'} = { params.environment_type};
    parameter_table{i,'tf'} = { mat2str(params.tf)};
    
end

    parameter_table{n_tot,'log10c0'} = {num2str(params.log10c0,16)};
    parameter_table{n_tot,'Omega'} ={ mat2str(params.Omega,16)};
    parameter_table{n_tot,'deltaE'} = { mat2str(params.deltaE,16)};
    parameter_table{n_tot,'delta_Omega2'} = { mat2str(params.delta_Omega2)};
    parameter_table{n_tot,'max_batches'} = { mat2str(params.max_batches)};
    parameter_table{n_tot,'environment_type'} = { 'deterministic'};
    parameter_table{n_tot,'tf'} = { mat2str(params.tf)};

parameter_table{:,'id'} = (1:n_tot)';

writetable(parameter_table,'automated_runs/parameters/check_stochastic_coexistence.csv','Delimiter',',')


%% Look at bernoulli stochastic coexistence


clear;clc
load('automated_runs/parameters/default_params_struct.mat');
table_columns = {'log10c0','deltaE','Omega','delta_Omega2','max_batches','environment_type','tf',};
n_par = length(table_columns);
n_tot = 1000;
params.deltaE = [0;0.05];
params.log10c0 = 2;
params.tf = 10*params.tf;
params.Omega = [0;1e-6 + (params.deltaE(2)./params.E(1)).*log(1+(10.^params.log10c0))/(params.rho0.*params.tf)];
params.delta_Omega2 = 0.5*params.Omega(2);
params.environment_type = 'stochastic-Omega-bernoulli';
params.max_batches = 2e5;

parameter_table = cell(n_tot,n_par);

parameter_table = cell2table(parameter_table,'VariableNames',table_columns);

for i = 1:(n_tot-1)
    
    parameter_table{i,'log10c0'} = {num2str(params.log10c0,16)};
    parameter_table{i,'Omega'} ={ mat2str(params.Omega,16)};
    parameter_table{i,'deltaE'} = { mat2str(params.deltaE,16)};
    parameter_table{i,'delta_Omega2'} = { mat2str(params.delta_Omega2)};
    parameter_table{i,'max_batches'} = { mat2str(params.max_batches)};
    parameter_table{i,'environment_type'} = { params.environment_type};
    parameter_table{i,'tf'} = { mat2str(params.tf)};
    
end

    parameter_table{n_tot,'log10c0'} = {num2str(params.log10c0,16)};
    parameter_table{n_tot,'Omega'} ={ mat2str(params.Omega,16)};
    parameter_table{n_tot,'deltaE'} = { mat2str(params.deltaE,16)};
    parameter_table{n_tot,'delta_Omega2'} = { mat2str(params.delta_Omega2)};
    parameter_table{n_tot,'max_batches'} = { mat2str(params.max_batches)};
    parameter_table{n_tot,'environment_type'} = { 'deterministic'};
    parameter_table{n_tot,'tf'} = { mat2str(params.tf)};

parameter_table{:,'id'} = (1:n_tot)';

writetable(parameter_table,'automated_runs/parameters/bernoulli_stochastic_coexistence.csv','Delimiter',',')

%% Look at critical stochastic coexistence

clear;clc
load('automated_runs/parameters/default_params_struct.mat');
table_columns = {'log10c0','deltaE','Omega','delta_Omega2','max_batches','environment_type','tf',};
n_par = length(table_columns);
n_tot = 1000;
params.deltaE = [0;0.05];
params.log10c0 = 2;
params.tf = 10*params.tf;
params.Omega = [0;(params.deltaE(2)./params.E(1)).*log(1+(10.^params.log10c0))/(params.rho0.*params.tf)];
params.delta_Omega2 = 0.5*params.Omega(2);
params.environment_type = 'stochastic-Omega';
params.max_batches = 2e5;

parameter_table = cell(n_tot,n_par);

parameter_table = cell2table(parameter_table,'VariableNames',table_columns);

for i = 1:(n_tot-1)
    
    parameter_table{i,'log10c0'} = {num2str(params.log10c0,16)};
    parameter_table{i,'Omega'} ={ mat2str(params.Omega,16)};
    parameter_table{i,'deltaE'} = { mat2str(params.deltaE,16)};
    parameter_table{i,'delta_Omega2'} = { mat2str(params.delta_Omega2)};
    parameter_table{i,'max_batches'} = { mat2str(params.max_batches)};
    parameter_table{i,'environment_type'} = { params.environment_type};
    parameter_table{i,'tf'} = { mat2str(params.tf)};
    
end

    parameter_table{n_tot,'log10c0'} = {num2str(params.log10c0,16)};
    parameter_table{n_tot,'Omega'} ={ mat2str(params.Omega,16)};
    parameter_table{n_tot,'deltaE'} = { mat2str(params.deltaE,16)};
    parameter_table{n_tot,'delta_Omega2'} = { mat2str(params.delta_Omega2)};
    parameter_table{n_tot,'max_batches'} = { mat2str(params.max_batches)};
    parameter_table{n_tot,'environment_type'} = { 'deterministic'};
    parameter_table{n_tot,'tf'} = { mat2str(params.tf)};

parameter_table{:,'id'} = (1:n_tot)';

writetable(parameter_table,'automated_runs/parameters/critical_stochastic_coexistence.csv','Delimiter',',')


%% Looking at high omega stochastic coexistence (10k)

clear;clc
load('automated_runs/parameters/default_params_struct.mat');
table_columns = {'log10c0','deltaE','Omega','delta_Omega2','max_batches','environment_type','tf',};
n_par = length(table_columns);
n_tot = 10000;
params.deltaE = [0;0.05];
params.log10c0 = 2;
params.tf = 10*params.tf;
params.Omega = [0;1e-6 + (params.deltaE(2)./params.E(1)).*log(1+(10.^params.log10c0))/(params.rho0.*params.tf)];
params.delta_Omega2 = 0.5*params.Omega(2);
params.environment_type = 'stochastic-Omega';
params.max_batches = 4e5;

parameter_table = cell(n_tot,n_par);

parameter_table = cell2table(parameter_table,'VariableNames',table_columns);

for i = 1:(n_tot-1)
    
    parameter_table{i,'log10c0'} = {num2str(params.log10c0,16)};
    parameter_table{i,'Omega'} ={ mat2str(params.Omega,16)};
    parameter_table{i,'deltaE'} = { mat2str(params.deltaE,16)};
    parameter_table{i,'delta_Omega2'} = { mat2str(params.delta_Omega2)};
    parameter_table{i,'max_batches'} = { mat2str(params.max_batches)};
    parameter_table{i,'environment_type'} = { params.environment_type};
    parameter_table{i,'tf'} = { mat2str(params.tf)};
    
end

    parameter_table{n_tot,'log10c0'} = {num2str(params.log10c0,16)};
    parameter_table{n_tot,'Omega'} ={ mat2str(params.Omega,16)};
    parameter_table{n_tot,'deltaE'} = { mat2str(params.deltaE,16)};
    parameter_table{n_tot,'delta_Omega2'} = { mat2str(params.delta_Omega2)};
    parameter_table{n_tot,'max_batches'} = { mat2str(params.max_batches)};
    parameter_table{n_tot,'environment_type'} = { 'deterministic'};
    parameter_table{n_tot,'tf'} = { mat2str(params.tf)};

parameter_table{:,'id'} = (1:n_tot)';

writetable(parameter_table,'automated_runs/parameters/high_Omega_coexistence_10K.csv','Delimiter',',')

%% Looking at weak Omega noise, but tf~9 not ~90 (10k)
% To be matched with the relevant tf noise
clear;clc
load('automated_runs/parameters/default_params_struct.mat');
table_columns = {'log10c0','deltaE','Omega','delta_Omega2','max_batches','environment_type','tf'};
n_par = length(table_columns);
n_tot = 10000;
params.deltaE = [0;0.05];
params.log10c0 = 2;
% params.tf = 10*params.tf;
params.Omega = [0;1e-5 + (params.deltaE(2)./params.E(1)).*log(1+(10.^params.log10c0))/(params.rho0.*params.tf)];
params.delta_Omega2 = 0.1*params.Omega(2);
params.environment_type = 'stochastic-Omega';
params.max_batches = 3e5;

parameter_table = cell(n_tot,n_par);

parameter_table = cell2table(parameter_table,'VariableNames',table_columns);

for i = 1:(n_tot-1)
    
    parameter_table{i,'log10c0'} = {num2str(params.log10c0,16)};
    parameter_table{i,'Omega'} ={ mat2str(params.Omega,16)};
    parameter_table{i,'deltaE'} = { mat2str(params.deltaE,16)};
    parameter_table{i,'delta_Omega2'} = { mat2str(params.delta_Omega2)};
    parameter_table{i,'max_batches'} = { mat2str(params.max_batches)};
    parameter_table{i,'environment_type'} = { params.environment_type};
    parameter_table{i,'tf'} = { mat2str(params.tf)};
    
end

    parameter_table{n_tot,'log10c0'} = {num2str(params.log10c0,16)};
    parameter_table{n_tot,'Omega'} ={ mat2str(params.Omega,16)};
    parameter_table{n_tot,'deltaE'} = { mat2str(params.deltaE,16)};
    parameter_table{n_tot,'delta_Omega2'} = { mat2str(params.delta_Omega2)};
    parameter_table{n_tot,'max_batches'} = { mat2str(params.max_batches)};
    parameter_table{n_tot,'environment_type'} = { 'deterministic'};
    parameter_table{n_tot,'tf'} = { mat2str(params.tf)};

parameter_table{:,'id'} = (1:n_tot)';

writetable(parameter_table,'automated_runs/parameters/small_Omega_coexistence_tf_notlong_correctedDrift_10K.csv','Delimiter',',')

%% Looking at weak Omega noise, but tf~4.5 not (half of ~9). Omega drift scaled accordingly (10k)
% To be matched with the relevant tf noise
clear;clc
load('automated_runs/parameters/default_params_struct.mat');
table_columns = {'log10c0','deltaE','Omega','delta_Omega2','max_batches','environment_type','tf'};
n_par = length(table_columns);
n_tot = 10000;
params.deltaE = [0;0.05];
params.log10c0 = 2;
params.tf = params.tf/2;
params.Omega = [0;2e-5 + (params.deltaE(2)./params.E(1)).*log(1+(10.^params.log10c0))/(params.rho0.*params.tf)];
params.delta_Omega2 = 0.1*params.Omega(2);
params.environment_type = 'stochastic-Omega';
params.max_batches = 3e5;

parameter_table = cell(n_tot,n_par);

parameter_table = cell2table(parameter_table,'VariableNames',table_columns);

for i = 1:(n_tot-1)
    
    parameter_table{i,'log10c0'} = {num2str(params.log10c0,16)};
    parameter_table{i,'Omega'} ={ mat2str(params.Omega,16)};
    parameter_table{i,'deltaE'} = { mat2str(params.deltaE,16)};
    parameter_table{i,'delta_Omega2'} = { mat2str(params.delta_Omega2)};
    parameter_table{i,'max_batches'} = { mat2str(params.max_batches)};
    parameter_table{i,'environment_type'} = { params.environment_type};
    parameter_table{i,'tf'} = { mat2str(params.tf)};
    
end

    parameter_table{n_tot,'log10c0'} = {num2str(params.log10c0,16)};
    parameter_table{n_tot,'Omega'} ={ mat2str(params.Omega,16)};
    parameter_table{n_tot,'deltaE'} = { mat2str(params.deltaE,16)};
    parameter_table{n_tot,'delta_Omega2'} = { mat2str(params.delta_Omega2)};
    parameter_table{n_tot,'max_batches'} = { mat2str(params.max_batches)};
    parameter_table{n_tot,'environment_type'} = { 'deterministic'};
    parameter_table{n_tot,'tf'} = { mat2str(params.tf)};

parameter_table{:,'id'} = (1:n_tot)';

writetable(parameter_table,'automated_runs/parameters/small_Omega_coexistence_tf_half_correctedDrift_10K.csv','Delimiter',',')


%% Looking at weak tf noise, but tf~9 not ~90 (10k)
% To be matched with the relevant tf noise
clear;clc
load('automated_runs/parameters/default_params_struct.mat');
table_columns = {'log10c0','deltaE','Omega','delta_tf','max_batches','environment_type','tf'};
n_par = length(table_columns);
n_tot = 10000;
params.deltaE = [0;0.05];
params.log10c0 = 2;
% params.tf = 10*params.tf;
params.Omega = [0;1e-5 + (params.deltaE(2)./params.E(1)).*log(1+(10.^params.log10c0))/(params.rho0.*params.tf)];
params.delta_tf = 0.1*params.tf;
params.environment_type = 'stochastic-tf';
params.max_batches = 3e5;

parameter_table = cell(n_tot,n_par);

parameter_table = cell2table(parameter_table,'VariableNames',table_columns);

for i = 1:(n_tot-1)
    
    parameter_table{i,'log10c0'} = {num2str(params.log10c0,16)};
    parameter_table{i,'Omega'} ={ mat2str(params.Omega,16)};
    parameter_table{i,'deltaE'} = { mat2str(params.deltaE,16)};
    parameter_table{i,'delta_tf'} = { mat2str(params.delta_tf)};
    parameter_table{i,'max_batches'} = { mat2str(params.max_batches)};
    parameter_table{i,'environment_type'} = { params.environment_type};
    parameter_table{i,'tf'} = { mat2str(params.tf)};
    
end

    parameter_table{n_tot,'log10c0'} = {num2str(params.log10c0,16)};
    parameter_table{n_tot,'Omega'} ={ mat2str(params.Omega,16)};
    parameter_table{n_tot,'deltaE'} = { mat2str(params.deltaE,16)};
    parameter_table{n_tot,'delta_tf'} = { mat2str(params.delta_tf)};
    parameter_table{n_tot,'max_batches'} = { mat2str(params.max_batches)};
    parameter_table{n_tot,'environment_type'} = { 'deterministic'};
    parameter_table{n_tot,'tf'} = { mat2str(params.tf)};

parameter_table{:,'id'} = (1:n_tot)';

writetable(parameter_table,'automated_runs/parameters/small_Omega_coexistence_stochastic_tf_notlong_correctedDrift_10K.csv','Delimiter',',')



%% Checking for opposite stochastic coexistence

clear;clc
load('automated_runs/parameters/default_params_struct.mat');
table_columns = {'log10c0','deltaE','Omega','delta_Omega2','max_batches','environment_type','tf',};
n_par = length(table_columns);
n_tot = 2000;
params.deltaE = [0;0.05];
params.log10c0 = 2;
params.tf = 10*params.tf;
params.Omega = [0;-1e-6 + (params.deltaE(2)./params.E(1)).*log(1+(10.^params.log10c0))/(params.rho0.*params.tf)];
params.delta_Omega2 = 0.5*params.Omega(2);
params.environment_type = 'stochastic-Omega';
params.max_batches = 2e5;

parameter_table = cell(n_tot,n_par);

parameter_table = cell2table(parameter_table,'VariableNames',table_columns);

for i = 1:(n_tot-1)
    
    parameter_table{i,'log10c0'} = {num2str(params.log10c0,16)};
    parameter_table{i,'Omega'} ={ mat2str(params.Omega,16)};
    parameter_table{i,'deltaE'} = { mat2str(params.deltaE,16)};
    parameter_table{i,'delta_Omega2'} = { mat2str(params.delta_Omega2)};
    parameter_table{i,'max_batches'} = { mat2str(params.max_batches)};
    parameter_table{i,'environment_type'} = { params.environment_type};
    parameter_table{i,'tf'} = { mat2str(params.tf)};
    
end

    parameter_table{n_tot,'log10c0'} = {num2str(params.log10c0,16)};
    parameter_table{n_tot,'Omega'} ={ mat2str(params.Omega,16)};
    parameter_table{n_tot,'deltaE'} = { mat2str(params.deltaE,16)};
    parameter_table{n_tot,'delta_Omega2'} = { mat2str(params.delta_Omega2)};
    parameter_table{n_tot,'max_batches'} = { mat2str(params.max_batches)};
    parameter_table{n_tot,'environment_type'} = { 'deterministic'};
    parameter_table{n_tot,'tf'} = { mat2str(params.tf)};

parameter_table{:,'id'} = (1:n_tot)';

writetable(parameter_table,'automated_runs/parameters/opposite_stochastic_coexistence.csv','Delimiter',',')

%% Deterministic and stochastic IC sweep

clear;clc
load('automated_runs/parameters/default_params_struct.mat');
table_columns = {'log10c0','deltaE','Omega','max_batches','environment_type','tf','rho_at_t0','delta_Omega2'};
n_par = length(table_columns);
n_tot = 1200;
params.deltaE = [0;0.05];
params.log10c0 = 2;
params.tf = 10*params.tf;
log_ratio_range = linspace(-30,30,200);
rho1_at_t0_range = 1./(1 + exp(log_ratio_range));
params.Omega = [0;1e-6 + (params.deltaE(2)./params.E(1)).*log(1+(10.^params.log10c0))/(params.rho0.*params.tf)];
params.delta_Omega2 = 0.5*params.Omega(2);
environment_type_cell = {'stochastic-Omega','stochastic-Omega',...
    'stochastic-Omega','stochastic-Omega',...
    'stochastic-Omega','full-run-deterministic'};
params.max_batches = 3e5;

parameter_table = cell(n_tot,n_par);

parameter_table = cell2table(parameter_table,'VariableNames',table_columns);
tn = 1;
for j = 1:length(environment_type_cell)
    for i = 1:length(log_ratio_range)
        parameter_table{tn,'log10c0'} = {num2str(params.log10c0,16)};
        parameter_table{tn,'Omega'} ={ mat2str(params.Omega,16)};
        parameter_table{tn,'rho_at_t0'} = {mat2str([rho1_at_t0_range(i);1-rho1_at_t0_range(i)])};
        parameter_table{tn,'deltaE'} = { mat2str(params.deltaE,16)};
        parameter_table{tn,'max_batches'} = { mat2str(params.max_batches)};
        parameter_table{tn,'environment_type'} = { environment_type_cell{j}};
        parameter_table{tn,'delta_Omega2'} = { mat2str(params.delta_Omega2)};
        parameter_table{tn,'tf'} = { mat2str(params.tf)};
        tn = tn + 1;
    end
end


parameter_table{:,'id'} = (1:n_tot)';

writetable(parameter_table,'automated_runs/parameters/deterministic_stochastic_IC_sweep.csv','Delimiter',',')

%% Checking for moderate omega timescales

clear;clc
load('automated_runs/parameters/default_params_struct.mat');
table_columns = {'log10c0','deltaE','Omega','delta_Omega2','max_batches','environment_type','tf',};
n_par = length(table_columns);
n_tot = 1000;
params.deltaE = [0;0.05];
params.log10c0 = 2;
params.tf = 10*params.tf;
params.Omega = [0;0.5e-6 + (params.deltaE(2)./params.E(1)).*log(1+(10.^params.log10c0))/(params.rho0.*params.tf)];
params.delta_Omega2 = 0.5*params.Omega(2);
params.environment_type = 'stochastic-Omega';
params.max_batches = 3e5;

parameter_table = cell(n_tot,n_par);

parameter_table = cell2table(parameter_table,'VariableNames',table_columns);

for i = 1:(n_tot-1)
    
    parameter_table{i,'log10c0'} = {num2str(params.log10c0,16)};
    parameter_table{i,'Omega'} ={ mat2str(params.Omega,16)};
    parameter_table{i,'deltaE'} = { mat2str(params.deltaE,16)};
    parameter_table{i,'delta_Omega2'} = { mat2str(params.delta_Omega2)};
    parameter_table{i,'max_batches'} = { mat2str(params.max_batches)};
    parameter_table{i,'environment_type'} = { params.environment_type};
    parameter_table{i,'tf'} = { mat2str(params.tf)};
    
end

    parameter_table{n_tot,'log10c0'} = {num2str(params.log10c0,16)};
    parameter_table{n_tot,'Omega'} ={ mat2str(params.Omega,16)};
    parameter_table{n_tot,'deltaE'} = { mat2str(params.deltaE,16)};
    parameter_table{n_tot,'delta_Omega2'} = { mat2str(params.delta_Omega2)};
    parameter_table{n_tot,'max_batches'} = { mat2str(params.max_batches)};
    parameter_table{n_tot,'environment_type'} = { 'full-run-deterministic'};
    parameter_table{n_tot,'tf'} = { mat2str(params.tf)};

parameter_table{:,'id'} = (1:n_tot)';

writetable(parameter_table,'automated_runs/parameters/moderate_Omega_stochastic.csv','Delimiter',',')


%% Looking at high omega stochastic coexistence with half noise (10k)

clear;clc
load('automated_runs/parameters/default_params_struct.mat');
table_columns = {'log10c0','deltaE','Omega','delta_Omega2','max_batches','environment_type','tf',};
n_par = length(table_columns);
n_tot = 10000;
params.deltaE = [0;0.05];
params.log10c0 = 2;
params.tf = 10*params.tf;
params.Omega = [0;1e-6 + (params.deltaE(2)./params.E(1)).*log(1+(10.^params.log10c0))/(params.rho0.*params.tf)];
params.delta_Omega2 = 0.1*params.Omega(2);
params.environment_type = 'stochastic-Omega';
params.max_batches = 4e5;

parameter_table = cell(n_tot,n_par);

parameter_table = cell2table(parameter_table,'VariableNames',table_columns);

for i = 1:(n_tot-1)
    
    parameter_table{i,'log10c0'} = {num2str(params.log10c0,16)};
    parameter_table{i,'Omega'} ={ mat2str(params.Omega,16)};
    parameter_table{i,'deltaE'} = { mat2str(params.deltaE,16)};
    parameter_table{i,'delta_Omega2'} = { mat2str(params.delta_Omega2)};
    parameter_table{i,'max_batches'} = { mat2str(params.max_batches)};
    parameter_table{i,'environment_type'} = { params.environment_type};
    parameter_table{i,'tf'} = { mat2str(params.tf)};
    
end

    parameter_table{n_tot,'log10c0'} = {num2str(params.log10c0,16)};
    parameter_table{n_tot,'Omega'} ={ mat2str(params.Omega,16)};
    parameter_table{n_tot,'deltaE'} = { mat2str(params.deltaE,16)};
    parameter_table{n_tot,'delta_Omega2'} = { mat2str(params.delta_Omega2,16)};
    parameter_table{n_tot,'max_batches'} = { mat2str(params.max_batches)};
    parameter_table{n_tot,'environment_type'} = { 'full-run-deterministic'};
    parameter_table{n_tot,'tf'} = { mat2str(params.tf,16)};

parameter_table{:,'id'} = (1:n_tot)';

writetable(parameter_table,'automated_runs/parameters/high_Omega_fifth_noise_10K.csv','Delimiter',',')


%% Set up for production variable K sweep

clear;clc
rng(666)
n_ind1 = 41;
n_ind2 = 2;
n_ind3 = 1001;

table_columns = {'log10c0','K','deltaE'};
n_par = length(table_columns);
n_tot = n_ind1*n_ind2*n_ind3;

load('automated_runs/parameters/default_params_struct.mat');

parameter_table = cell(n_tot,n_par);

parameter_table = cell2table(parameter_table,'VariableNames',table_columns);

log10c0_range = linspace(-3,3,n_ind1);
deltaE2_range = [0.05,1];

Nr_den_sq = compute_Nr_den_sq(log10c0_range,params.rho0,1,0);

tn = 1;
for i = 1:n_ind1
    log10c0 = log10c0_range(i);
    for k = 1:n_ind2
        deltaE2 = deltaE2_range(k);
        
        predicted_deltaK = ...
            deltaE2_range(k).*log((10.^log10c0_range(i) + params.rho0)./params.rho0)...
            ./(Nr_den_sq(i)./params.E(1));

        K2_range = linspace(1+0.1*predicted_deltaK,1+6*predicted_deltaK,n_ind3);
        
        for j = 1:n_ind3
            parameter_table{tn,'log10c0'} = {num2str(log10c0,16)};
            parameter_table{tn,'K'} ={ mat2str([1;K2_range(j)],16)};
            parameter_table{tn,'deltaE'} = { mat2str([0;deltaE2],16)};
            
            tn = tn + 1;
        end
    end
    disp(i)
end

parameter_table{:,'id'} = (1:n_tot)';

parameter_table = parameter_table(randperm(n_tot),:);

writetable(parameter_table,'automated_runs/parameters/production_variable_K_sweep.csv','Delimiter',',')


%% Checking for mean deviation with tf noise

clear;clc
load('automated_runs/parameters/default_params_struct.mat');
table_columns = {'log10c0','deltaE','Omega','delta_tf','max_batches','environment_type','tf'};
n_par = length(table_columns);
n_tot = 3000;
params.deltaE = [0;0.05];
params.log10c0 = 2;
params.tf = params.tf;
params.Omega = [0;1e-6 + (params.deltaE(2)./params.E(1)).*log(1+(10.^params.log10c0)/params.rho0)/(params.tf)];
params.delta_tf = 0.5*params.tf;
params.environment_type = 'stochastic-tf';
params.max_batches = 3e5;

parameter_table = cell(n_tot,n_par);

parameter_table = cell2table(parameter_table,'VariableNames',table_columns);

for i = 1:(n_tot-1)
    
    parameter_table{i,'log10c0'} = {num2str(params.log10c0,16)};
    parameter_table{i,'Omega'} ={ mat2str(params.Omega,16)};
    parameter_table{i,'deltaE'} = { mat2str(params.deltaE,16)};
    parameter_table{i,'delta_tf'} = { mat2str(params.delta_tf,16)};
    parameter_table{i,'max_batches'} = { mat2str(params.max_batches)};
    parameter_table{i,'environment_type'} = { params.environment_type};
    parameter_table{i,'tf'} = { mat2str(params.tf,16)};
    
end

    parameter_table{n_tot,'log10c0'} = {num2str(params.log10c0,16)};
    parameter_table{n_tot,'Omega'} ={ mat2str(params.Omega,16)};
    parameter_table{n_tot,'deltaE'} = { mat2str(params.deltaE,16)};
    parameter_table{n_tot,'delta_tf'} = { mat2str(params.delta_tf,16)};
    parameter_table{n_tot,'max_batches'} = { mat2str(params.max_batches)};
    parameter_table{n_tot,'environment_type'} = { 'full-run-deterministic'};
    parameter_table{n_tot,'tf'} = { mat2str(params.tf,16)};

parameter_table{:,'id'} = (1:n_tot)';

writetable(parameter_table,'automated_runs/parameters/tf_noise_deviation_check.csv','Delimiter',',')




%% Checking for tf noise with skewed initial abundances

clear;clc
load('automated_runs/parameters/default_params_struct.mat');
table_columns = {'log10c0','deltaE','rho_at_t0','Omega','delta_tf','max_batches','environment_type','tf'};
n_par = length(table_columns);
n_tot = 3000;
params.deltaE = [0;0.05];
params.log10c0 = 2;
params.rho_at_t0 = [1-1e-5;1e-5];
params.tf = params.tf;
params.Omega = [0;1e-6 + (params.deltaE(2)./params.E(1)).*log(1+(10.^params.log10c0)/params.rho0)/(params.tf)];
params.delta_tf = 0.5*params.tf;
params.environment_type = 'stochastic-tf';
params.max_batches = 3e5;

parameter_table = cell(n_tot,n_par);

parameter_table = cell2table(parameter_table,'VariableNames',table_columns);

for i = 1:(n_tot-1)
    
    parameter_table{i,'log10c0'} = {num2str(params.log10c0,16)};
    parameter_table{i,'Omega'} ={ mat2str(params.Omega,16)};
    parameter_table{i,'deltaE'} = { mat2str(params.deltaE,16)};
    parameter_table{i,'rho_at_t0'} = { mat2str(params.rho_at_t0,16) }; 
    parameter_table{i,'delta_tf'} = { mat2str(params.delta_tf,16)};
    parameter_table{i,'max_batches'} = { mat2str(params.max_batches)};
    parameter_table{i,'environment_type'} = { params.environment_type};
    parameter_table{i,'tf'} = { mat2str(params.tf,16)};
    
end

    parameter_table{n_tot,'log10c0'} = {num2str(params.log10c0,16)};
    parameter_table{n_tot,'Omega'} ={ mat2str(params.Omega,16)};
    parameter_table{n_tot,'deltaE'} = { mat2str(params.deltaE,16)};
    parameter_table{n_tot,'rho_at_t0'} = { mat2str(params.rho_at_t0,16) }; 
    parameter_table{n_tot,'delta_tf'} = { mat2str(params.delta_tf,16)};
    parameter_table{n_tot,'max_batches'} = { mat2str(params.max_batches)};
    parameter_table{n_tot,'environment_type'} = { 'full-run-deterministic'};
    parameter_table{n_tot,'tf'} = { mat2str(params.tf,16)};

parameter_table{:,'id'} = (1:n_tot)';

writetable(parameter_table,'automated_runs/parameters/tf_noise_different_starting.csv','Delimiter',',')