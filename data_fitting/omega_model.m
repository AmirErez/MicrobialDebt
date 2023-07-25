function [model_survival_frac] = omega_model(t_vec,non_debt_init,wd,wn) 

%This script models the death kinetics of mixed debtor-non-debtor
%population

model_survival_frac = (1-non_debt_init).*exp(-wd*t_vec) + non_debt_init.*exp(-wn*t_vec);

end