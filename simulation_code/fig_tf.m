%% 

E = 1;
log10c0_vec = linspace(-4,4,60);
gtype_cell = {'linear','hill','hill'};
h_vec = [NaN,1,2];


for i = 1:length(log10c0_vec)
    for j = 1:length(h_vec)
        [tf(i,j),Qf(i,j)] = ...
            compute_single_species_tf(0.99,log10c0_vec(i),1,1,gtype_cell{j},h_vec(j));
    end
end

newfigure(3,2);
hold on
hmap=['#00429d'; '#4c76b2'; '#a5ab58'];
plot(10.^log10c0_vec,tf(:,1),'--', 'Color', hmap(1,:), 'LineWidth',1.5,'DisplayName','Linear')
plot(10.^log10c0_vec,tf(:,2),'-', 'Color', hmap(2,:), 'LineWidth',1.5,'DisplayName','Monod')
plot(10.^log10c0_vec,tf(:,3),'--', 'Color', hmap(3,:), 'LineWidth',1.5,'DisplayName','Hill-2')
set(gca,'XScale','log','YScale','linear')
xlabel('$c_0/K$', 'Interpreter','latex')
ylabel('99\% consumption: $E\, t^*$', 'Interpreter','latex')
ylim([1,10])
xlim(10.^[-4,4]);
xticks(10.^[-4,0,4]);
tf = compute_single_species_tf(0.99,4,1,1,'hill',1);
plot([10.^log10c0_vec(1), 10.^log10c0_vec(end)], [tf, tf],':k')
text(10^3, 8.5, '$t_f$', 'Interpreter','latex')
l=legend('Interpreter','latex','Location','northwest');
l.String=l.String(1:3);
l.Position = [0.1841    0.5609    0.3171    0.2597];
print(gcf, '-dpng','../figures/fig_tf.png','-r600');
print(gcf, '-dsvg','../figures/fig_tf.svg');
