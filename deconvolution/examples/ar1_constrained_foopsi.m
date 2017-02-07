%% test modulate for all oasis functions. 
col = {[0 114 178],[0 158 115], [213 94 0],[230 159 0],...
    [86 180 233], [204 121 167], [64 224 208], [240 228 66]}; % colors
plot_cvx = false; 

%% example 7:  constrained-foopsi, AR1
g = 0.95;         % AR coefficient 
noise = .3; 
T = 3000; 
framerate = 30;     
firerate = 0.5; 
b = 0;              % baseline 
N = 1;              % number of trials 
seed = 13;          % seed for genrating random variables 
[y, true_c, true_s] = gen_data(g, noise, T, framerate, firerate, b, N, seed); 

% cvx solution 
[c_cvx, s_cvx] = constrained_foopsi(y, g, noise); 
% case 1: all parameters are known 
[c_oasis, s_oasis] = deconvolveCa(y, 'ar1', g, 'constrained', 'sn', noise);  %#ok<*ASGLU>

figure('name', 'constrained-FOOPSI, AR1, known: g, sn', 'papersize', [15, 4]);
plot_cvx = true; 
show_results; 
plot_cvx = false; 

% case 2: nothing is known, estimate g with auto-correlation method
[c_oasis, s_oasis,options] = deconvolveCa(y, 'ar1', 'constrained'); 

fprintf('true gamma:        %.3f\n', g); 
fprintf('estimated gamma:   %.3f\n', options.pars); 

figure('name', 'FOOPSI, AR1, estimated: g, sn', 'papersize', [15, 4]); 
show_results; 

% case 3: nothing is know, estimate g with auto-correlation method first
% and then update it to minimize the RSS
[c_oasis, s_oasis, options] = deconvolveCa(y, 'ar1', 'constrained', ...
    'optimize_pars'); 

fprintf('true gamma:        %.3f\n', g); 
fprintf('estimated gamma:   %.3f\n', options.pars); 

figure('name', 'FOOPSI, AR1, estimated: g, sn, update:g', 'papersize', [15, 4]); 
show_results; 

% case 4: nothing is know, estimate g with auto-correlation method first
% and then update it to minimize the RSS, the baseline is also unknown
true_b = 0.5; 
[c_oasis, s_oasis, options] = deconvolveCa(y+true_b, 'ar1', g,...
    'constrained','optimize_b', 'sn', noise); 
fprintf('true gamma:        %.3f\n', g); 
fprintf('estimated gamma:   %.3f\n', options.pars); 
fprintf('true b:       %.3f\n', true_b); 
fprintf('estimated b:       %.3f\n', options.b); 
fprintf('tuning parameter:  %.3f\n', options.lambda); 

figure('name', 'FOOPSI, AR1, estimated: g, sn, lambda', 'papersize', [15, 4]); 
show_results; 

% case 5: nothing is know, estimate g with auto-correlation method first
% and then update it to minimize the RSS, the baseline is also unknown
true_b = 0.5; 
[c_oasis, s_oasis, options] = deconvolveCa(y+true_b, 'ar1',...
    'constrained','optimize_b', 'optimize_pars'); 
fprintf('true gamma:        %.3f\n', g); 
fprintf('estimated gamma:   %.3f\n', options.pars); 
fprintf('estimated b:       %.3f\n', options.b); 
fprintf('tuning parameter:  %.3f\n', options.lambda); 

figure('name', 'FOOPSI, AR1, estimated: g, sn, lambda, update:g', 'papersize', [15, 4]); 
show_results; 

%% 