function tau_dr = ar2exp(g)
%% get parameters of the convolution kernel for AR2 process 
temp = roots([1, -g(1), -g(2)]);
d = max(temp); 
r = min(temp);
tau_d = -1/log(d); 
tau_r = -1/log(r); 

tau_dr = [tau_d, tau_r]; 