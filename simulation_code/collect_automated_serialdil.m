
function collect_automated_serialdil(param_file, sim_name,num_jobs,run_per_job)

%Import parameter table
parameter_table = readtable(param_file);

%Loop through files and fill storage structures
for i = 1:num_jobs

    %Load file if it exists
    sim_file = ['automated_runs/results/',sim_name,'/sim_group_',num2str(i),'.mat'];
    if isfile(sim_file)
        load(sim_file)
        if exist('output_cell')
            %Add in data
            for j = 1:run_per_job
                output = output_cell{j};
                true_ind = run_per_job*(i-1) + j;
                table_ind = parameter_table.id == true_ind;
                parameter_table.final_rho{table_ind} = {mat2str(output.rho(:,end))};
                parameter_table.n_batches{table_ind} = {mat2str(size(output.rho,2))};
                parameter_table.end_S{table_ind} = {mat2str(output.ShannonS(end))};
                parameter_table.end_NutIntegrals{table_ind} = {mat2str(output.NutIntegrals{end})};
            end
            disp(['Collected group ',num2str(i)])

        else
            for j = 1:run_per_job
                %Put in nans if unfinished
                true_ind = run_per_job*(i-1) + j;
                table_ind = parameter_table.id == true_ind;
                parameter_table.final_rho{table_ind} = {mat2str([NaN;NaN])};
                parameter_table.n_batches{table_ind} = {mat2str(NaN)};
                parameter_table.end_S{table_ind} = {mat2str(NaN)};
                parameter_table.end_NutIntegrals{table_ind} = {mat2str(NaN)};
            end
        end
    else
        for j = 1:run_per_job
            %Put in nans if unfinished
            true_ind = run_per_job*(i-1) + j;
            table_ind = parameter_table.id == true_ind;
            parameter_table.final_rho{table_ind} = {mat2str([NaN;NaN])};
            parameter_table.n_batches{table_ind} = {mat2str(NaN)};
            parameter_table.end_S{table_ind} = {mat2str(NaN)};
            parameter_table.end_NutIntegrals{table_ind} = {mat2str(NaN)};
        end
    end
end

writetable(parameter_table,['automated_runs/results/',sim_name,'_results_table.csv'])

disp(['Saved collected results to ',['automated_runs/results/',sim_name,'_results_table.csv']])

end
