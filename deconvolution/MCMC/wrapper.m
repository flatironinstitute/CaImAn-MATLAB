clear;
load traceData.mat;
addpath utilities
addpath(genpath('../constrained-foopsi'));
%Y = squeeze(traceData.traces(129,7,:));   % pick a particular trace (low SNR)
Y = mean(squeeze(traceData.traces(:,7,:))); % average over ROI (high SNR)

%% run MCMC sampler and plot results
params.p = 1;
params.print_flag = 1;
params.B = 300;
SAMP = cont_ca_sampler(Y,params);
plot_continuous_samples(SAMP,Y);