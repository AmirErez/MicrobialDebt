function average_stochastic_runs(sim_name,num_jobs,run_per_job)

%Loop through files and fill storage structures
tn = 1;
for i = 1:num_jobs
    
    %Load file if it exists
    sim_file = ['results/',sim_name,'/sim_group_',num2str(i),'.mat'];
    if isfile(sim_file)
        load(sim_file)
        if exist('output_cell')
            %Add in data
            for j = 1:run_per_job
		disp(tn)
                output = output_cell{j};
                if ~exist('running_sum')
                    running_sum = output.rho;
                else
                    running_sum = running_sum + output.rho;
                end
                final_rho{tn} = output.rho(:,end);
                tn = tn + 1;
            end
        end
    end
end

average_rho = running_sum/(tn-1);

save(['results/',sim_name,'_averaged_results.mat'],'final_rho','average_rho')

end

