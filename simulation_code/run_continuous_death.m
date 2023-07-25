%% Run new numerical code

clear;clc

params.p = 1;
params.m = 2;
params.K = 1;
params.max_batches = 1e6;
params.E = [1;1];
params.P = 1;
params.d = [0;0];
params.rho0 = 1;
params.rho_at_t0 = [1/2;1/2];
params.psi = [1;1.01];
params.Omega = [1;1];
params.alpha = [1;1];
params.errtype = 1;
params.dilution_method = 'floating';
plt = 0;

log10c0_vec = -1;

tic
for j = 1:length(log10c0_vec)
    params.log10c0 = log10c0_vec(j);
    %params.batch_length = 2*log(10.^log10c0_vec(j));
    %params.batch_length = 10*log(10.^params.log10c0);
    params.batch_length = 20;
    %predicted_omega = (10.^params.log10c0).*(params.psi(2)-1)/params.batch_length;
    predicted_omega = log(10.^params.log10c0).*(params.psi(2)-1)/params.batch_length;
    omega2_vec = linspace(0.1*predicted_omega,10*predicted_omega,40);


    for i = 1:length(omega2_vec)
        params.omega = [0;omega2_vec(i)];
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
save('local_sim_data/continuous_check.mat')


%%
H = @(x) nansum(-x.*log2(x));

for j = 1:length(log10c0_vec) 
    for i = 1:length(omega2_vec)
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
plot(omega2_vec,H_mat(:,i),'.-','Color',colors(i,:),...
    'DisplayName',num2str(log10c0_vec(i)),'LineWidth',1);
%'DisplayName',num2str(log2((10.^c0_vec(i) + params.rho0)./params.rho0)));
end
%xlim([1,50])
%ylim([1e-7,1])
legend()
ylabel('H')
xlabel('\omega')

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


