function output = serialdil_odesolver(params, plt)

% Amir Erez 2022-02-08.
% This function simulates serial dilutions
% Simulation ends dependent on extinction tolerance and RAE tolerance

output = struct;
output.rho = zeros(params.m, params.max_batches);
output.NutIntegrals = cell(params.max_batches,1);
tol = 1e-7;
er=Inf;
params.rho_initial = params.rho_at_t0;
params.mean_log10c0 = params.log10c0;
params.mean_tf = params.tf;


temp_log10c0s = [2,1,2.2];

for cnt=1:params.max_batches  
    if(mod(cnt,1000)==0)
       disp(['   Done ' num2str(cnt) '  ;  err ' num2str(max(er))]);
       disp(['Populations are ',mat2str(params.rho_at_t0./sum(params.rho_at_t0))])
    end
    
    if strcmp(params.environment_type,'stochastic-c0')
        logn_mode_c0 = log(10.^params.mean_log10c0) - (params.log10c0_var^2)/2;
%         params.log10c0 = log10(lognrnd(logn_mode_c0,params.log10c0_var));
        params.log10c0 = temp_log10c0s(cnt);
    elseif strcmp(params.environment_type,'stochastic-tf')
        %logn_mode_tf = log(params.mean_tf) - (params.tf_var^2)/2;
        %params.tf = lognrnd(logn_mode_tf,params.tf_var);
        params.tf = 2*params.mean_tf*rand(1);
    elseif strcmp(params.environment_type,'stochastic-both')
        logn_mode_c0 = log(10.^params.mean_log10c0) - (params.log10c0_var^2)/2;
        params.log10c0 = log10(lognrnd(logn_mode_c0,params.log10c0_var));
        %logn_mode_tf = log(params.mean_tf) - (params.tf_var^2)/2;
        %params.tf = lognrnd(logn_mode_tf,params.tf_var);
        params.tf = 2*params.mean_tf*rand(1);
    end
    
    %Simulate batch
    [rho_sigma, Nr, c_i, t] = multispeciesbatch_odesolver(params,0);
    
    rho_sigma(rho_sigma<1e-20) = 0;
    
    %Choose dilution method and perform
    params.rho_at_t0 = rho_sigma(:,end).*exp(-params.Omega.*params.tf).*...
        ((params.rho0)./(params.rho0 + 10.^params.mean_log10c0));
    output.tot_rho_post_catast(cnt) = sum(rho_sigma(:,end).*exp(-params.Omega.*params.tf));
       
    output.rho(:,cnt) = params.rho_at_t0;
    output.NutIntegrals{cnt} = Nr';
    output.Q(:,cnt) = (10.^params.log10c0 - c_i(:,end))./(10.^params.log10c0); 
    output.t(cnt) = params.tf;
    output.log10c0(cnt) = params.log10c0;
    
    if sum(params.rho_at_t0) == 0
        disp('All species extinct, ending simulation')
        break
    end
    
    %Compute relative error between batches of populations or diversity
    if(cnt>1) 
        if params.errtype == 1
            er = abs((output.rho(:,cnt) - output.rho(:,cnt-1)))/sum(output.rho(:,cnt));
        else
            Dcurr = -sum(output.rho(:,cnt).*log(output.rho(:,cnt)));
            Dprev = -sum(output.rho(:,cnt-1).*log(output.rho(:,cnt-1)));
            er = abs((Dcurr - Dprev)/Dcurr);
        end
    
        if (max(er) < tol) && ~contains(params.environment_type,'stochastic') %Exit loop of error below threshold
            break;
%             disp(['The model has been ended due to an RAE of '...
%                    num2str(transpose(er))])
%             disp('The populations are')
%             disp(bstore(:,i+1))
        end
%         if min(bstore(:,i+1)) < extol %Exit loop if an organism has died
%             o = 1;
%             disp(['The model has been ended due to an population of '...
%                   num2str(min(bstore(:,i+1)))])
%         end
    end
end

if(cnt==params.max_batches) && ~contains(params.environment_type,'stochastic')
    disp('Reached max batches. Increase params.max_batches to attain steady state.')
end

%PROCESS DATA--------------------------------------------------------------

output.rho = output.rho(:,1:cnt);
output.NutIntegrals = output.NutIntegrals(1:cnt);
output.ShannonS = calc_entropy_nats(output.rho);

%PLOTTING------------------------------------------------------------------

if plt == 1
    
    %PLOT END BIOMASS RATIO DYNAMICS
    colors = jet(params.m);
    figure
    set(gcf, 'Position', [500 250 750 600]);
    for h=1:params.m  %Plot end-batch biomass ratios vs. transfer number
        hold on
        if params.p ==3
        semilogy(output.rho(h,:), 'LineWidth',1.5,'color',params.alpha(h,:))
        else 
        semilogy(output.rho(h,:), 'LineWidth',1.5,'color',colors(h,:))
        end
    end
    %title(['Competition between ' num2str(params.m) ' species for ' num2str(params.p) ...
    %    ' nutrients at ' num2str(params.P(1)) '/' num2str(params.P(2))])
    xlabel('Transfer number')
    ylabel('Population fraction at batch start')
    leg = [];
    for h = 1:params.m %Generate legend from strategies
        legi = [];
        for k = 1:params.p
            if k < params.p
                legi = [legi num2str(params.alpha(h,k),'%5.2f') '/'];
            else
                legi = [legi num2str(params.alpha(h,k),'%5.2f')];
            end
        end
        leg = [leg; legi];
    end
    legend(leg)
    set(gca,'YScale', 'log')
    ylim([1e-10,100])
    %SIMPLEX PLOTS IF P = 3
    if params.p == 3 %Plots this only when there are three nutrients
        figure %Plot strategies and nutrient supply on the simplex
        set(gcf, 'Position', [500 500 750 600]);
        for i = 1:params.m
        ternplot(params.alpha(i,1),params.alpha(i,2),params.alpha(i,3),'o','MarkerEdgeColor',params.alpha(i,:),'majors',...
            0,'MarkerFaceColor',params.alpha(i,:),'majors',0)
        hold on
        end
        ternplot(params.P(1),params.P(2),params.P(3),'kd','majors',0,'LineWidth',2)
    end
    
    
end

end

%© 2020 GitHub, Inc.
