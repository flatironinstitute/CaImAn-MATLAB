%% function 
clear; clc; close all; 
T = 1000; 
s = double(rand(1, T)>0.98); 
sig = 0.04; 

% example
tau_d = 5; 
tau_r = 1; 
nMax = 100; 
pars = [tau_d, tau_r]; 
kernel = create_kernel('exp2', pars, nMax); 

t = 1:kernel.nMax; 
gt = kernel.fhandle(kernel.pars, t);  
c1 = conv(s, gt); 
c1 = c1(1:T); 
y1 = c1 + randn(1, T) * sig; 

%% use the true convolution kernel  
kernel0 = kernel; 
kernel0.pars = [10, 0.1]; 
kernel0.bound_pars = false; 
figure('position', [1,1,1500, 200]); 
plot(y1); 
hold on; 
plot(c1, 'r', 'linewidth', 2); 
plot(-s*0.1, 'r', 'linewidth', 2); 
tic; 
[chat, shat, kernel_fit, iters] = deconvCa(y1, kernel0, 2, true, false);
toc; 
plot(chat,'-.g','linewidth', 2); 
plot(-shat*0.1, '-.g', 'linewidth', 2); %, '-.'); 
legend('data', 'ground truth: c','ground truth: s', 'OASIS:c', 'OASIS: s'); % 'gound truth: s', 'OASIS: c', 'OASIS: s');