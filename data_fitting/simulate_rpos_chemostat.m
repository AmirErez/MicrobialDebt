function [x_nd_vec,tout_cell,yout_cell] = simulate_rpos_chemostat(target_tau_vec,E,DeltaE_01,DeltaE_03,DeltaE_06,debtor_init,Gamma,delta_vec)

%Right-hand side function of ODE
% RHS = @(t,y) [y(1)*E*y(3) - delta*y(1); ...
%     (E+DeltaE)*y(3)*y(2) - delta*y(2);...
%     Gamma - y(1)*E*y(3) - (E+DeltaE)*y(3)*y(2) - delta*y(3)];


%Get the target time from the target tau
target_t_vec = target_tau_vec./delta_vec*log(2);

for i = 1:length(target_tau_vec)
  
    %Establish RHS function
    if length(delta_vec) == 1
        delta = delta_vec;
    else
        delta = delta_vec(i);
    end

    if delta == 0.1
        DeltaE = DeltaE_01;
    elseif delta == 0.3
        DeltaE = DeltaE_03;
    elseif delta == 0.6
        DeltaE = DeltaE_06;
    else
        disp('Improper delta specified.')
        DeltaE = 1;
    end

    RHS = @(t,y) [y(1)*E*Gamma/(y(1)*E + y(2)*(E+DeltaE)) - delta*y(1); ...
       y(2)*(E+DeltaE)*Gamma/(y(1)*E + y(2)*(E+DeltaE)) - delta*y(2)];

    %Set up initial conditions
    y0 = [Gamma/delta debtor_init/(1-debtor_init)*Gamma/delta];

    %Simulate
    [tout_cell{i},yout_cell{i}] = ode45(RHS,[0,target_t_vec(i)],y0);

    %Compute relative abundance
    x_nd_vec(i,1) = yout_cell{i}(end,1)/sum(yout_cell{i}(end,:));

end

x_nd_vec = log10(x_nd_vec);

end