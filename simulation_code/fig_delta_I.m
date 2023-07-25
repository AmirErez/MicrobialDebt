% Plots the relation between the normal Monod integral
% and the one with (c+K)^2 in the denominator. Used to explain
% The variable K coexistance.

log10c0s = -3:0.1:3;
K=1;
r0=1;
c0s = 10.^log10c0s;
dilutions = r0./(r0+c0s);

[Nr_den_sq] = compute_Nr_den_sq(log10c0s,r0,K,0);


newfigure(3,2);

loglog(c0s, Nr_den_sq,'-k', 'LineWidth',1.5, 'DisplayName', '$\Delta_I / \Delta_K$')
%semilogx(c0s, 1./exp(Nr_den_sq),'-k', 'LineWidth',2,'DisplayName','$\exp\left[\int_0^{t_f}\frac{c}{(c+K)^2}\right]$')

legend('Interpreter','latex', 'Location','northwest')
hold on

loglog(c0s, log(1+c0s/r0),'--r', 'LineWidth',1.5, 'DisplayName','$c_0\ll K$')
% semilogx(c0s, 1./(1+c0s/r0),'--r', 'LineWidth',2, 'DisplayName','$c_0\ll K$')


%vals = ( -(c0s+K+r0).*(c0s./(c0s+K)) + (c0s+r0).*log((c0s+K).*(c0s+r0)/K/r0) ) ./ ((c0s+K+r0).^2);
%vals = ( -(c0s+K+r0).*(c0s./(c0s+K)) + (c0s+r0).*log((c0s+K).*(c0s+r0)/K/r0) ) ./ ((c0s+K+r0).^2);
% vals = (log(c0s.^2/K/r0)-1)./c0s;
vals = (log(c0s.^2/K/r0))./c0s;

% loglog(c0s, vals, '--b', 'LineWidth',2)
% semilogx(c0s, exp(-vals), '--b', 'LineWidth',2, 'DisplayName', '$c_0\gg K$')
% semilogx(c0s, c0s.^(-2./c0s), '--g', 'LineWidth',2, 'DisplayName', '$c_0\gg K$')
%semilogx(c0s, (1./(1+c0s/r0)).^(2./c0s),'--g', 'LineWidth',2, 'DisplayName','$c_0\gg K$')
plot(c0s, log(((1+c0s/r0)).^(2./c0s)),'--b', 'LineWidth',1.5, 'DisplayName','$c_0\gg K$')


ylim([0,1]); 
xlim(10.^[-3,3])
xticks(10.^[-3,0,3]);
xlabel('$c_0/K$', 'Interpreter','latex');
ylabel('$\Delta_I / \Delta_K$', 'Interpreter','latex');
set(gca,'YScale','linear')
print(gcf,'-dpng', ['../figures/response_dk_delta_I.png'], '-r600')
print(gcf,'-dsvg', ['../figures/response_dk_delta_I.svg'])