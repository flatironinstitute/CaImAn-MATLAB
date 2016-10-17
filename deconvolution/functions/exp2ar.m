function g = exp2ar(tau_dr)
%% convert a convolution function to an AR(2) model 

d = exp(-1/tau_dr(1)); 
r = exp(-1/tau_dr(2));

g(1) = d+r; 
g(2) = -d*r; 