%% test modulate for all oasis functions. 
col = {[0 114 178],[0 158 115], [213 94 0],[230 159 0],...
    [86 180 233], [204 121 167], [64 224 208], [240 228 66]}; % colors
plot_cvx = false; 


%% threshold foopsi, convolution kernel  
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
temp = roots([1, -g(1), -g(2)]);
d = max(temp); 
r = min(temp);
w = 200;
ht = (exp(log(d)*(1:w)) - exp(log(r)*(1:w))) / (d-r); % convolution kernel
  
% case 1: all parameters are known, the kernel is the sum of two
% exponential functions
smin = 0.5; 
pars = [d, r]; 
[c_oasis, s_oasis] = deconvolveCa(y, 'exp2', pars, 'thresholded', 'smin', smin);  %#ok<*ASGLU>
figure('name', 'threshold, exp2, known: taur, taud, smin', 'papersize', [15, 4]); 
show_results; 

% case 1: all parameters are known 
smin = 0.5; 
[c_oasis, s_oasis] = deconvolveCa(y, 'kernel', ht, 'thresholded', 'smin', smin);  %#ok<*ASGLU>
figure('name', 'threshold, kernel, known: kernel, smin', 'papersize', [15, 4]); 
show_results; 
%%%%%%%%%%%%%%  END %%%%%%%%%%%%%%%%%%
