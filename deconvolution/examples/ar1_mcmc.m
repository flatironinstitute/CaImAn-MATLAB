%% test modulate for all oasis functions. 
col = {[0 114 178],[0 158 115], [213 94 0],[230 159 0],...
    [86 180 233], [204 121 167], [64 224 208], [240 228 66]}; % colors
plot_cvx = false; 

%% example 1:  foopsi, AR1 model. This model is used when the sampling rate is low
g = 0.95;         % AR coefficient 
noise = .3; 
T = 3000; 
framerate = 30;     
firerate = 0.5; 
b = 0;              % baseline 
N = 1;              % number of trials 
seed = 13;          % seed for genrating random variables 
[y, true_c, true_s] = gen_data(g, noise, T, framerate, firerate, b, N, seed); 

% case 1: all parameters are known 
lambda = 2.4; 
params.p = 1; 
params.B = 1; 
[c_oasis, s_oasis] = deconvolveCa(y, 'mcmc', params);  %#ok<*ASGLU>
[c_cvx, s_cvx] = foopsi(y, g, lambda); 

figure('name', 'FOOPSI, AR1, known: g, lambda', 'papersize', [15, 4]); 
plot_cvx = true; 
show_results; 
plot_cvx = false; 
