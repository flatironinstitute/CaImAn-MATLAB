function ht = exp2kernel(taus, nmax)
%% create the convolution kernel given tau_rise and tau_decay 

%% inputs: 
%   taus:   1*2 vector, [tau_decay tau_rise]
%   nmax:   scalar, length of the kernel 

%% outputs: 
%   ht:     nmax*1 kernel, convolution kernel 

%% Author: Pengcheng Zhou, Carnegie Mellon University, 2016 

t = (1:nmax)'; 
d = exp(-1./taus(1)); 
r = exp(-1./taus(2)); 
ht = (exp(log(d)*t) - exp(log(r)*t) ) / (d - r); 
ht = ht/sum(ht); 