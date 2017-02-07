function tau = tau_d2c(g,dt)

% convert discrete time constantswith resolution dt to continuous 
% h(t) = (1-exp(-t/tau(1)))*exp(-t/tau(2))

% Author: Eftychios A. Pnevmatikakis, 2016, Simons Foundation

gr = max(roots([1,-g(:)']),0);
p1_continuous = log(min(gr))/dt; 
p2_continuous = log(max(gr))/dt;
tau_1 = -1/p1_continuous;                   %tau h - smaller (tau_d * tau_r)/(tau_d + tau_r)
tau_2 = -1/p2_continuous;                   %tau decay - larger

tau_rise = 1/(1/tau_1 - 1/tau_2);
tau = [tau_rise,tau_2];                     %tau_h , tau_d