%% test modulate for all oasis functions. 
col = {[0 114 178],[0 158 115], [213 94 0],[230 159 0],...
    [86 180 233], [204 121 167], [64 224 208], [240 228 66]}; % colors
plot_cvx = false; 


%% threshold, AR2 model 
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
% case 1: all parameters are known 
smin = 0.5; 
[c_oasis, s_oasis] = deconvolveCa(y, 'ar2', g, 'thresholded', 'smin', smin);  %#ok<*ASGLU>
figure('name', 'threshold, AR2, known: g, smin', 'papersize', [15, 4]); 
show_results; 

% case 2: know smin
smin = 0.5; 
[c_oasis, s_oasis, options] = deconvolveCa(y, 'ar2', 'sn', noise, 'thresholded',...
    'smin', smin); 
fprintf('true gamma:        %.3f\t %.3f\n', g(1), g(2)); 
fprintf('estimated gamma:   %.3f\t %.3f\n', options.pars(1),  options.pars(2)); 

figure('name', 'threshold, AR2, known:smin, estimated: g', 'papersize', [15, 4]); 
show_results; 

%% case 3: know smin, update g
smin = 0.5; 
[c_oasis, s_oasis, options] = deconvolveCa(y, 'ar2', 'sn', noise, 'thresholded',...
    'smin', smin, 'optimize_pars'); 
fprintf('true gamma:        %.3f\t %.3f\n', g(1), g(2)); 
fprintf('estimated gamma:   %.3f\t %.3f\n', options.pars(1),  options.pars(2)); 

figure('name', 'threshold, AR2, known:smin, estimated: g', 'papersize', [15, 4]); 
show_results; 

%% case 3: estimate smin 
[c_oasis, s_oasis, options] = deconvolveCa(y, 'ar2', 'sn', noise, 'thresholded',...
    'optimize_smin','optimize_pars', 'thresh_factor', 1); 
% fprintf('true gamma:        %.3f\t %.3f\n', g(1), g(2)); 
% fprintf('estimated gamma:   %.3f\t %.3f\n', options.pars(1),  options.pars(2)); 
fprintf('estimated smin:    %.3f\n', options.smin); 
figure('name', 'threshold, AR2, known:smin, estimated: g', 'papersize', [15, 4]); 
show_results; 
%%%%%%%%%%%%%%  END %%%%%%%%%%%%%%%%%%
