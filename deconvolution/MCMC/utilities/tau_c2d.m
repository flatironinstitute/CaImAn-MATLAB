function [g,h1] = tau_c2d(tau_r,tau_d,dt)

% convert continuous time constants to discrete with resolution dt
% h(t) = (1-exp(-t/tau_r))*exp(-t/tau_d)
% g: discrete time constants
% h1: h(dt);
% h*s can be written in discrete form as filter(h1,[1,-g],s)

% Author: Eftychios A. Pnevmatikakis, 2016, Simons Foundation

A = [-(2/tau_d+1/tau_r), - (tau_r+tau_d)/(tau_r*tau_d^2); 1 0];
lc = eig(A*dt);
ld = exp(lc);
g = [sum(ld),-prod(ld)];
h1 = (1-exp(-dt/tau_r))*exp(-dt/tau_d);