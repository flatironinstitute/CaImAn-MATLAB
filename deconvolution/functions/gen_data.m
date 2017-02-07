function [Y, truth, trueSpikes] = gen_data(gam, noise, T, framerate, ...
    firerate, b, N, seed)

%% input arguments 
if ~exist('gam', 'var') || isempty(gam)
    gam = .95; 
end 
if ~exist('noise', 'var') || isempty(noise)
    noise = .3; 
end 
if ~exist('T', 'var') || isempty(T)
    T = 3000; 
end 
if ~exist('framerate', 'var') || isempty(framerate)
    framerate = 30; 
end 
if ~exist('firerate', 'var') || isempty(firerate)
    firerate = .5; 
end 
if ~exist('b', 'var') || isempty(b)
    b = 0; 
end 
if ~exist('N', 'var') || isempty(N)
    N = 20; 
end 
if ~exist('seed', 'var') || isempty(seed)
    seed = 13;  
end 

%% run simulation 
rng(seed); 
trueSpikes = (rand(N, T) < firerate/framerate); 
truth = double(trueSpikes); 
p = length(gam); 
gam = [flipud(reshape(gam, [], 1)); 1]; 

for t=(p+1):T
    truth(:, t) = truth(:, (t-p):t) * gam;
end 

Y = b + truth + noise * randn(N, T); 
        




















