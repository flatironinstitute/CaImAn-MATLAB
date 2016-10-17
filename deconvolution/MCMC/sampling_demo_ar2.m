% creating a demo for the MCMC sampler

clearvars;
addpath utilities
dt = 1e-1;   % Data is generated with provided dt but sampled with Dt = 1
T = 7000;
ld = 0.05;   % rate spikes per second

s = rand(1,round(T/dt)) < ld*dt;
tau_rise = 2;
tau_decay = 10;
hmax = tau_decay/(tau_decay+tau_rise)*(tau_rise/(tau_decay+tau_rise))^(tau_rise/tau_decay);
[g,h1] = tau_c2d(tau_rise,tau_decay,dt);

b = hmax/4;                                                         % baseline value
cin = [.2*b,.15*b];                                                 % initial concentration
c = [cin,filter(h1,[1,-g],s(3:end),filtic(h1,[1,-g],cin))] + b;     % true calcium at original resolution
c_true = c(round(1/dt):round(1/dt):round(T/dt));                    % true calcium at sampled resolution
sg = hmax/4;                                                        % noise level
y = c_true + sg*randn(1,length(c_true));                            % noisy observations

figure;subplot(2,1,1); plot(dt:dt:T,c); hold all; scatter(1:T,y,'r*'); 
                       legend('True Calcium','Observed Values');

%%  Run constrained deconvolution approach (constrained foopsi)

[g2,h2] = tau_c2d(tau_rise,tau_decay,1);   % generate discrete time constants

[ca_foopsi,cb,c1,~,~,spikes_foopsi] = constrained_foopsi(y,[],[],g2);

% do some plotting
spiketimes{1} =  find(s)*dt;
spikeRaster = samples_cell2mat(spiketimes,T,1);
f = find(spikeRaster);
spikes = zeros(sum(spikeRaster),2);
count = 0;
for cnt = 1:length(f)
    spikes(count+(1:spikeRaster(f(cnt))),1) = f(cnt);
    spikes(count+(1:spikeRaster(f(cnt))),2) = 0.95 + 0.025*max(spikes_foopsi)*(1:spikeRaster(f(cnt)));
    count = count + spikeRaster(f(cnt));
end

subplot(2,1,2);
stem(spikes_foopsi); hold all; 
    scatter(spikes(:,1),spikes(:,2)-0.95+max(spikes_foopsi),15,'magenta','filled');
    axis([1,T,0,max(spikes(:,2))-0.95+max(spikes_foopsi)]);
    title('Foopsi Spikes','FontWeight','bold','Fontsize',14); xlabel('Timestep','FontWeight','bold','Fontsize',16);
    legend('Foopsi Spikes','Ground Truth');
    drawnow;

%% run continuous time MCMC sampler

params.p = 2;
params.g = g2;
params.sp = spikes_foopsi;   % pass results of foopsi for initialization (if not, they are computed)
params.c = ca_foopsi;
params.b = cb;
params.c1 = c1;
params.sn = sg;
params.marg = 0;

SAMPLES = cont_ca_sampler(y,params);    %% MCMC   
plot_continuous_samples(SAMPLES,y(:));
M = plot_marginals(SAMPLES.ss,T,spikeRaster);