function spikeRaster = samples_cell2mat(sampleCell,T,Dt)

% Constructs matrix of spike raster plot from cell array of spike times

% Inputs:
% sampleCell:   Cell array with spike times in continuous time (SAMPLES.ss)
% T:            End time of trace/experiment
% Dt:           Time-bin resolution (default: 1)

% Output:
% spikeRaster:  Spike raster plot matrix

% Author: Eftychios A. Pnevmatikakis and Josh Merel

if nargin == 2
    Dt = 1;
end
bins = 0:Dt:(T-Dt);
nsamples = length(sampleCell);
spikeRaster = zeros(nsamples,length(bins));
for i = 1:nsamples
    tmp = histc([sampleCell{i}(:); inf],[bins, (T+1)]);
    spikeRaster(i,:) = tmp(1:(end-1))';
end