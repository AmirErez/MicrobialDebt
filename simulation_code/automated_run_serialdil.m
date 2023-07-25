
function automated_run_serialdil(default_params_file,table_file,group_num,results_folder,run_per_job)

rng('shuffle')
if(isstr(group_num)), group_num=str2num(group_num);end
if(isstr(run_per_job)),run_per_job=str2num(run_per_job);end
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
%Generate results file and save
results_file = [results_folder,'/','sim_group_',num2str(group_num),'.mat'];
precollected_results_file = [results_folder,'/','precollected_sim_group_',num2str(group_num),'.mat'];

if(exist(results_file))
   disp(['File ' results_file ' exists. Not overwriting, quitting'])
   quit(0);
end
plt = 0;

%Loop through points
for i = start_end(1):start_end(2)
    %Only run simulations that are required
    if i <= size(params_table,1)
        %Read in default parameter structure
        load(default_params_file);
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

        true_ind = i - start_end(1) + 1;

        disp(['Running simulation ',num2str(i),'.'])

        %Run simulation and put in cell
        output_cell{true_ind} = serialdil_odesolver(params,plt);
        params_cell{true_ind} = params;

        %Make precollected data
        batches = 1:100:size(output_cell{true_ind}.rho,2);
        precollected_cell{true_ind}.batches = batches;
        precollected_cell{true_ind}.rho1 = output_cell{true_ind}.rho(1,batches);
        precollected_cell{true_ind}.rho2 = output_cell{true_ind}.rho(2,batches);
    end

end

save(results_file,'output_cell','params_cell');
save(precollected_results_file,'precollected_cell', 'params_cell');

disp(['Raw results saved to ',results_file,'.'])


end
