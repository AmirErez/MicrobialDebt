clear;clc
mkdir -p ../figures

%This figure generates figures 1A and 1C

%Load in data from Notley-McRobb
fig1_data = readtable('Notley-McRobb_data.xlsx','Sheet','Figure 1 Notley-McRobb 2002');
fig1_data = fig1_data(strcmp(fig1_data.limitingNutrient,'glucose'),:);
fig1_data.x_OfRpos_(fig1_data.x_OfRpos_ > 100) = 100;
fig1_data.x_debtor = fig1_data.x_OfRpos_/100;
xData = fig1_data.time_generations_;
yData = log10(fig1_data.x_debtor);
delta_vec = fig1_data.dilutionRate_h__1_;

%Define fit
ft = fittype('simulate_rpos_chemostat(target_tau_vec,E,DeltaE_01,DeltaE_03,DeltaE_06,debtor_init,Gamma,delta_vec)',...
    'coefficients',{'E','DeltaE_01','DeltaE_03','DeltaE_06','debtor_init','Gamma'},'independent',{'target_tau_vec'},...
    'problem',{'delta_vec'});
opts = fitoptions( 'Method', 'NonlinearLeastSquares');
opts.Lower = [1,0,0,0,0,1];
opts.StartPoint = [1 0.1 0.1 0.1 1e-2 1];
n_upper = 1e3;
opts.Upper = [1,n_upper,n_upper,n_upper,1,1];

%Fit model to data
[fitresult, gof] = fit( xData, yData, ft, opts,'problem',{delta_vec});
c_ints = confint(fitresult);
[x_nd_vec,tout_cell,yout_cell] = simulate_rpos_chemostat(xData,fitresult.E,...
    fitresult.DeltaE_01,fitresult.DeltaE_03,fitresult.DeltaE_06,fitresult.debtor_init,fitresult.Gamma,delta_vec);

%Display fitresult
fitresult


%Run to get confidence intervals
uplow = predint(fitresult,xData,0.95,'functional','on');
uplow = 10.^uplow;
uplow(uplow < 0) = 0;
fill_between_lines = @(X1,X2,Y1,Y2,C) fill( [transpose(X1) fliplr(transpose(X2))],  [transpose(Y1) fliplr(transpose(Y2))], C,'EdgeColor','none','FaceAlpha',0.2);

%% Plot fig 1A (evolution trajectories)

colors = [230 159 0; 0 114 178; 0 158 115]./255;
unique_delta = unique(delta_vec);
newfigure(6,2);
hold on
FontSize = 11;
for i = 1:length(unique_delta)
    delta_ind = delta_vec == unique_delta(i);
    x = xData(delta_ind);
    y = 10.^yData(delta_ind);

    %Simulate equivalent trajectory
    theory_x = linspace(min(x),max(x),30);
    theory_delta = zeros(size(theory_x)) + unique_delta(i);
    x_nd_vec = simulate_rpos_chemostat(theory_x,fitresult.E,...
        fitresult.DeltaE_01,fitresult.DeltaE_03,fitresult.DeltaE_06,fitresult.debtor_init,fitresult.Gamma,theory_delta);
    theory_y = 10.^x_nd_vec;

    up_data = uplow(delta_ind,2);
    low_data = uplow(delta_ind,1);
    h = fill_between_lines(x,x, up_data,low_data,colors(i,:));
    set( get( get( h, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );

    fit_name = ['\delta = ',num2str(unique_delta(i))];
    h1(i) = plot(x,y,'o','MarkerFaceColor',colors(i,:),'Color',colors(i,:),'DisplayName',fit_name);
    h2(i) = plot(theory_x,theory_y,'-','LineWidth',1.5,'Color',colors(i,:),'HandleVisibility','off');
end
set(gca,'YScale','log')
xlabel('Generations','Interpreter','latex')
ylabel('Non-debtor fraction','Interpreter','latex')
legend(h1,'FontSize',FontSize)
ylim([1e-4,1])
set(gca,'FontSize',FontSize)

print(gcf,'../figures/fig1A.svg','-dsvg','-r600');


%% Plot fig1C (growth rate vs. deltaE)

newfigure(3,2);
hold on
FontSize = 11;
x = [0.1 0.3 0.6];
y = [fitresult.DeltaE_01, fitresult.DeltaE_03,fitresult.DeltaE_06];
err = y - c_ints(1,2:4);
plot(1./x,y,'ko','MarkerFaceColor','k');
errorbar(1./x,y,err,'k-','LineWidth',1.5);
set(gca,'XScale','linear','YScale','linear')
xlabel('Inverse dilution, $1/\delta$ (h)','Interpreter','latex')
ylabel({'Debtor mutant', 'advantage, $\Delta_E/E$'},'Interpreter','latex')
set(gca,'FontSize',FontSize)

print(gcf,'../figures/fig1C.svg','-dsvg','-r600');

