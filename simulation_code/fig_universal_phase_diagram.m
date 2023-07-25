clear;clc

n = 200;
log10c0_range = linspace(-3,3,n);
deltaE2 = 0.05;
E1 = 1;
rho0 = 1;
tf = 9.2018;
lower_y = 1e-5;
upper_y = 1e3;
newfigure(3,2);
debtor_color = [206 37 123]/255;
nondebtor_color = [15 104 194]/255;
hold on
c0_range = 10.^log10c0_range;
Omega = (deltaE2./E1).*log(1+c0_range./rho0)./tf;
chemostat_Omega = (deltaE2./E1).*c0_range./rho0./tf;

lower_bounds = zeros(size(c0_range)) + lower_y;
upper_bounds = zeros(size(c0_range)) + upper_y;
x = [c0_range,fliplr(c0_range)];
y1 = [lower_bounds,fliplr(Omega/deltaE2)];
y2 = [Omega/deltaE2,fliplr(upper_bounds)];
fill(x,y1,debtor_color,'LineWidth',1.5,'EdgeColor','k')
fill(x,y2,nondebtor_color,'LineWidth',1.5,'EdgeColor','k')
plot(c0_range,chemostat_Omega/(deltaE2),'--','LineWidth',1.5,'Color','k');
plot(c0_range,Omega/(deltaE2),'k-','LineWidth',1.52);
text(10.^(1.8),1e-3,'Debtor','Interpreter','latex',...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'Color','w')
text(10.^(-1.8),1e0,'Non-debtor','Interpreter','latex',...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'Color','w')
text(35,12,'Chemostat','Interpreter','latex',...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'Color','w','Rotation',26.5)

set(gca,'XScale','log','YScale','log')
set(gca,'TickDir','out')
set(gca,'YminorTick','off')
set(gca,'XminorTick','off')
xlabel('$c_0/K$','Interpreter','latex')
ylabel('$\Omega/\Delta_E$, or $\omega/\Delta_E$','Interpreter','latex')
xticks(10.^[-3,0,3]);
ylim([lower_y,upper_y])
yticks([1e-5,1e-1,1e3])
print(gcf, '-dpng','../figures/fig_universal_phase_diagram.png','-r600');
print(gcf, '-dsvg','../figures/fig_universal_phase_diagram.svg');
