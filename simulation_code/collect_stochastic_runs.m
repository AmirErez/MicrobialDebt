function collect_stochastic_runs(sim_name,num_jobs,run_per_job,b_range,b_interval)

%Loop through files and fill storage structures
tn = 1;
batches = b_range(1):b_interval:b_range(2);
for i = 1:num_jobs
    
    %Load file if it exists
    sim_file = ['automated_runs/results/',sim_name,'/sim_group_',num2str(i),'.mat'];
    if isfile(sim_file)
        try
           load(sim_file)
        catch
           disp(['Skipping ' sim_file]);
            continue;
        end
        if exist('output_cell')
            %Add in data
            for j = 1:run_per_job
                true_ind = run_per_job*(i-1) + j;
                disp([num2str(true_ind) ' ' sim_file])
                output = output_cell{j};
                rho1_mat(tn,:) = output.rho(1,batches);
                rho2_mat(tn,:) = output.rho(2,batches);
                true_ind_vec(tn) = true_ind;
                tn = tn + 1;
            end
        end
    end
end

save(['automated_runs/results/',sim_name,'_collected_results.mat'],'rho1_mat','rho2_mat','batches','true_ind_vec')

end

