function collect_precollected_stochastic_runs(sim_name,num_jobs,run_per_job,b_range,b_interval)

%Loop through files and fill storage structures
tn = 1;
batches = b_range(1):b_interval:b_range(2);
for i = 1:num_jobs
    
    %Load file if it exists
    sim_file = ['automated_runs/results/',sim_name,'/precollected_sim_group_',num2str(i),'.mat'];
    if isfile(sim_file)
        load(sim_file)
        if exist('precollected_cell')
            %Add in data
            for j = 1:run_per_job
                disp(tn)
                precollected_final.rho1_mat(tn,:) = precollected_cell{j}.rho1;
                precollected_final.rho2_mat(tn,:) = precollected_cell{j}.rho2;
                precollected_final.batch_cell{tn} = precollected_cell{j}.batches;
                tn = tn + 1;
            end
        end
    end
end

save(['automated_runs/results/',sim_name,'_precollected_results.mat'],'precollected_final')

end

