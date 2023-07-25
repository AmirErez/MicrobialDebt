%% 

E = 1;
log10c0_vec = linspace(-4,4,60);
c0_vec = 10.^log10c0_vec;
h_vec = [NaN,1,2];
newfigure(3,2);
set(gca, 'FontSize',11)
hold on

hmap=['#00429d'; '#4c76b2'; '#a5ab58'];
style={'--','-','--'};
leg={'Linear', 'Monod', 'Hill 2'};

for ii = 1:length(h_vec)
    if isnan(h_vec(ii))
        yvals = c0_vec;
    else
        yvals = (c0_vec.^h_vec(ii))./(c0_vec.^h_vec(ii)+1);
    end
    plot(10.^log10c0_vec, yvals, 'LineStyle',style{ii},'LineWidth',1.5,'Color',hmap(ii,:),...
        'DisplayName',leg{ii});
end


set(gca,'XScale','log','YScale','log')
xlabel('$c/K$', 'Interpreter','latex')
ylabel('Utilization, $g[c]$', 'Interpreter','latex')

%ylim([1,10])
xlim(10.^[-4,4]);
xticks(10.^[-4,0,4]);
text(10^(-3), 10^(-7), '$\sim \left(\frac{c}{K}\right)^2$', 'Interpreter','latex')
text(10^(1), 10^(3.8), '$\sim \left(\frac{c}{K}\right)$', 'Interpreter','latex')

l=legend('Interpreter','latex','Location','southeast');
l.String=l.String(1:3);
print(gcf, '-dpng','../figures/fig_growth_functions.png','-r600');
print(gcf, '-dsvg','../figures/fig_growth_functions.svg');
