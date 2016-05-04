function [srt,n_events,n_act] = order_components(YrA,C,sn,options)

% order components by how unlikely it is to produce the data under a no
% spiking hypothesis. The function counts how many times the trace is above
% the sum of the mode and a given number of standard deviations for a given
% number of consecutive intervals. 

% INPUTS

% YrA:      Noise for each temporal component
% C:        Temporal components
% sn:       Standard deviation of each component (optional, estimated if not given)
% options:  Options structure
%   options.nsd:    number of standard deviations to set threshold
%   options.nfr:    number of consecutive frames for considering an event

% OUPUTS

% srt:      Sorted components in decreasing order
% n_event:  Number of events for each components in decreasing order
% n_act:    Number of components with events (The srt(1:n_act) is the set of active components)

% Written by:
% Eftychios A. Pnevmatikakis, Simons Foundation, 2016

if nargin < 3 || isempty(sn)
    sn = std(YrA,[],2);
end

defoptions = CNMFSetParms;
if nargin < 4 || isempty(options)
    options = defoptions;
end

if ~isfield(options,'nsd') || isempty(options.nsd); options.nsd = defoptions.nsd; end  % number of sd above the mode
if ~isfield(options,'nfr') || isempty(options.nfr); options.nsd = defoptions.nfr; end  % number of consecutive frames

[K,T] = size(C);
CY = C + YrA;
md = zeros(K,1);
%bd = zeros(K,1);

for i = 1:K
    [bandwidth,density,xmesh]=kde(CY(i,:),2^floor(log2(T)-1));
    [~,id] = max(density);
    md(i) = xmesh(id);
    %bd(i) = bandwidth;
end

FF = false(size(CY));
FF2 = false(size(CY));
for i = 1:K
    ff = (CY(i,:)>md(i)+options.nsd*sn(i));
    FF(i,:) = ff;
    FF2(i,:) = bwareaopen(ff,options.nfr);
end

[n_events,srt] = sort(sum(FF2,2),'descend');
n_act = sum(n_events>0);