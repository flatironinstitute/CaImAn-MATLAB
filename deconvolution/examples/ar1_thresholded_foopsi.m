%% test modulate for all oasis functions. 
col = {[0 114 178],[0 158 115], [213 94 0],[230 159 0],...
    [86 180 233], [204 121 167], [64 224 208], [240 228 66]}; % colors
plot_cvx = false; 

%% example 4: hard threshold, AR1 model
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
smin = 0.5; 
[c_oasis, s_oasis, options] = deconvolveCa(y, 'ar1', g, 'thresholded', 'smin', smin);  %#ok<*ASGLU>

figure('name', 'threshold, AR1, known: g, lambda, smin', 'papersize', [15, 4]); 
show_results; 

% case 2: know smin
[c_oasis, s_oasis, options] = deconvolveCa(y, 'ar1', 'thresholded', 'smin', smin); 

fprintf('true gamma:        %.3f\n', g); 
fprintf('estimated gamma:   %.3f\n', options.pars); 

figure('name', 'threshold, AR1, known:smin, estimated: g', 'papersize', [15, 4]); 
show_results; 

% case 3: optimize the thershold, g, and the baseline
[c_oasis, s_oasis, options] = deconvolveCa(y, 'ar1', g,  ...
    'thresholded', 'optimize_smin', 'optimize_pars', 'thresh_factor', 0.99);  %#ok<*ASGLU>

figure('name', 'threshold, AR1, known: g, sn, estimate: smin, g', 'papersize', [15, 4]); 
show_results; 