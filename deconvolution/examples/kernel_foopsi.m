%% test modulate for all oasis functions. 
col = {[0 114 178],[0 158 115], [213 94 0],[230 159 0],...
    [86 180 233], [204 121 167], [64 224 208], [240 228 66]}; % colors
plot_cvx = false; 

%% example 3: foopsi, convolution kernel 
g = [1.7, -0.712];         % AR coefficient 
noise = 1; 
T = 3000; 
framerate = 30;     
firerate = 0.5; 
b = 0;              % baseline 
N = 20;              % number of trials 
seed = 3;          % seed for genrating random variables 
[Y, trueC, trueS] = gen_data(g, noise, T, framerate, firerate, b, N, seed); 
y = Y(1,:); 
true_c = trueC(1,:);  %#ok<*NASGU>
true_s = trueS(1,:);
taus = ar2exp(g); 
w = 200;
taus = ar2exp(g); 
ht = exp2kernel(taus, w); 
% case 1: use the difference of two exponential functions to construct a
% kernel 
lambda = 25; 
[c_oasis, s_oasis] = deconvolveCa(y, 'exp2', taus, 'foopsi', 'lambda', lambda, ...
    'shift', 100, 'window', 200);  %#ok<*ASGLU>

figure('name', 'FOOPSI, exp2, known: g, lambda', 'papersize', [15, 4]); 
show_results; 

% case 2: use the kernel directly 
lambda = 25; 
[c_oasis, s_oasis] = deconvolveCa(y, 'kernel', ht, 'foopsi', 'lambda', ...
    lambda, 'shift', 100, 'window', 200);  %#ok<*ASGLU>

figure('name', 'FOOPSI, kernel, known: g, lambda', 'papersize', [15, 4]); 
show_results; 


%% case 3: estimate the time constants 
lambda = 0; 
taus = ar2exp(g); 
[c_oasis, s_oasis, options] = deconvolveCa(y, 'exp2', 'foopsi', 'lambda', lambda, ...
    'shift', 100, 'window', 200, 'smin', 0.5);  %#ok<*ASGLU>

figure('name', 'FOOPSI, exp2, known: g, lambda', 'papersize', [15, 4]); 
show_results; 















