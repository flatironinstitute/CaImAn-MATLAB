function [M,fig] = plot_marginals(sampleCell,T,truespikes,int_show)

% Plots marginal posterior empirical pdfs for # of spikes for each timebin
% similar to figure 1B in Pnevmatikakis et al., Neuron 2016

% Inputs:
% sampleCell:   Cell array with spike times in continuous time (SAMPLES.ss)     
% T:            Number of timebins
% truespikes:   Number of true spikes per timebin (vector of size T x 1)
% int_show:     Show only a specified interval (default: [1,T])

% Output:
% M:            matrix of empirical posterior pdfs for each timebin

% Author: Eftychios A. Pnevmatikakis, 2016, Simons Foundation

if nargin < 4
    int_show = 1:T;
end
nT = length(int_show);

if nargin == 2
    truespikes = -0.5*ones(1,nT);
end

if length(truespikes) == T
    truespikes = truespikes(int_show);
end

Mat = samples_cell2mat(sampleCell,T);
Mat = Mat(:,int_show);

mS = min(max([Mat(:);truespikes(:)])+1,6);
M = zeros(mS,nT);
for i = 1:nT
    M(:,i) = hist(Mat(:,i),0:mS-1)/size(Mat,1);
end

cmap = bone(100);
cmap(2:26,:) = [];
fig = figure;
imagesc(M); axis xy; 
%set(gca,'Ytick',[0.5:(mS+.5)],'Yticklabel',[-1:(mS)]);  %hold all; plot(traceData.spikeFrames+1); set(gca,'YLim',[1,5])
%set(gca,'YLim',[1.25,mS+.25]);
hold all; scatter(1:nT,truespikes+1,[],'m');
set(gca,'YLim',[1.5,mS+.125]);
pos = get(gca,'Position');
set(gca,'Ytick',0.5+[-0.5:mS-.5],'Yticklabel',[-1+(0:mS)]);
colormap(cmap);
ylabel('# of Spikes ','fontweight','bold','fontsize',14);
xlabel('Timestep ','fontweight','bold','fontsize',14);
title('Posterior Spike Histogram (MCMC) ','fontweight','bold','fontsize',14);
cbar = colorbar('Location','East');
cpos = get(cbar,'Position');
set(cbar,'Position',[pos(1)+pos(4),cpos(2:4)]);
%set(gca,'Xtick',[])
%set(cbar,'Color',[1,1,1]);
set(cbar,'Fontsize',12);