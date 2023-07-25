clear;clc

%This script generates figure 1B using data from figure 3 of Notley-McRobb
%2002

%Non-debtor initial abundance on days 1 2 3 4, from paper
rpos_abun = [1; .24; .15; .005]; 

%Load in data from Notley-McRobb
fig3_data = readtable('Notley-McRobb_data.xlsx','Sheet','Figure 3 Notley-McRobb 2002');
fig3_data.time_minutes_(fig3_data.time_minutes_ <0) = 0;
fig3_data.survivalPercentage(fig3_data.survivalPercentage <0) = 0;
fig3_data.survival_frac = fig3_data.survivalPercentage/100;
stressor = 'temp';
fig3_data = fig3_data(strcmp(fig3_data.stressor,stressor),:);
fig3_data.non_debt_init = rpos_abun(fig3_data.day);
t_vec = fig3_data.time_minutes_;
surv_vec = fig3_data.survival_frac;
non_debt_init = fig3_data.non_debt_init;

%Define fit
ft = fittype('omega_model(t_vec,non_debt_init,wd,wn)','coefficients',...
    {'wd','wn'},'independent',{'t_vec'},'problem',{'non_debt_init'});
opts = fitoptions( 'Method', 'NonlinearLeastSquares');
opts.Lower = [0,0];
opts.StartPoint = [0.1, 0.1];
opts.Upper = [100,100];


%Fit model to data
[fitresult, gof] = fit( t_vec, surv_vec, ft, opts,'problem',{non_debt_init});
c_ints = confint(fitresult);
confints = [fitresult.wd,fitresult.wn] -c_ints;
model_survival_frac = omega_model(t_vec,non_debt_init,fitresult.wd,fitresult.wn);


%Run to get confidence intervals
uplow = predint(fitresult,t_vec,0.95,'functional','on'); 
uplow(uplow < 0) = 0;


%Plot the fit and data
fill_between_lines = @(X1,X2,Y1,Y2,C) fill( [transpose(X1) fliplr(transpose(X2))],  [transpose(Y1) fliplr(transpose(Y2))], C,'EdgeColor','none','FaceAlpha',0.2);

newfigure(3,2);
hold on
colors = ['r','g','b','k'];
day_vec = unique(fig3_data.day);
for i = 1:length(day_vec)

    day_i_index = fig3_data.day == day_vec(i);

    t = t_vec(day_i_index);
    up_data = uplow(day_i_index,2);
    low_data = uplow(day_i_index,1); 
    h = fill_between_lines(t,t, up_data,low_data,colors(i));
    set( get( get( h, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    
    DisplayName = ['Day ',num2str(i), ' - $x_N = ',num2str(rpos_abun(i)),'$'];
    plot(t,surv_vec(day_i_index),'o','Color',colors(i),'MarkerFaceColor',colors(i),'DisplayName',DisplayName);
    h = plot(t,model_survival_frac(day_i_index),'-','Color',colors(i),'LineWidth',1.5,'DisplayName','');
    set( get( get( h, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
end

ylabel('Viable CFUs','Interpreter','latex')
xlabel('Time (minutes)','Interpreter','latex')
xlim([0,max(t_vec)+1])
h = legend('Interpreter','latex');
pos = get(h,'Position');
pos(4) = 0.1*pos(4);
set(h,'Position',pos)
set(gca,'FontSize',11)

print(gcf,'fig1B.svg','-dsvg','-r600');
