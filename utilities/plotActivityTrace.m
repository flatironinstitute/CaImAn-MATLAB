function plotActivityTrace( inferred, filtered, raw, plotControl)

% This code plots the activity trace for all the ROIs
% Author: Weijian Yang

% input parameters:
% inferred, filtered, raw (from the output of "signalExtraction.m")
% plotControl.frameRate: frame rate in fps
% plotControl.normalization: normalize the amplitude of the temporal trace for plotting purpose
% plotControl.sep: vertical (amplitude) offset between adjacent temporal traces 
% plotControl.displayLabel: whether to display ROI ID for the temporal traces 
% plotControl.plotInferred: whether the inferred signal is plotted
% plotControl.plotFiltered: whether the filtered signal is plotted
% plotControl.plotRaw: whether the raw signal is plotted
% plotControl.rollingView: whether the display is set such that only 20 temporal traces are displayed at a time, followed by any key to move onto the next 20 traces, and etc.

if isfield(plotControl,'frameRate') frameRate=plotControl.frameRate;
else frameRate=1; end

if isfield(plotControl,'normalization') normalization=plotControl.normalization;
else normalization=1; end

if isfield(plotControl,'sep') sep=plotControl.sep;
else sep=1.2; end

if isfield(plotControl,'displayLabel') displayLabel=plotControl.displayLabel;
else displayLabel=1; end

if isfield(plotControl,'plotInferred') plotInferred=plotControl.plotInferred;
else plotInferred=1; end

if isfield(plotControl,'plotFiltered') plotFiltered=plotControl.plotFiltered;
else plotFiltered=1; end

if isfield(plotControl,'plotRaw') plotRaw=plotControl.plotRaw;
else plotRaw=0; end

if isfield(plotControl,'rollingView') rollingView=plotControl.rollingView;
else rollingView=0; end

t=1/frameRate*(1:size(inferred,2));
K=size(inferred,1);

figure;
hold on;

if normalization==1
    for idx=1:K
        if plotRaw  plot(t,raw(idx,:)/max(raw(idx,:))-idx*sep,'color',[0.8 0.8 0.8],'linewidth',2); end
        if plotFiltered  plot(t,filtered(idx,:)/max(filtered(idx,:))-idx*sep,'color',[1 0.7 0.7],'linewidth',1.5); end
        if plotInferred  plot(t,inferred(idx,:)/max(inferred(idx,:))-idx*sep,'color',[0 0 1],'linewidth',1); end
    end
    ylabel('Normalized \DeltaF/F');
else
    for idx=1:K
        if plotRaw  plot(t,raw(idx,:)-idx*sep,'color',[0.8 0.8 0.8],'linewidth',2); end
        if plotFiltered  plot(t,filtered(idx,:)-idx*sep,'color',[1 0.7 0.7],'linewidth',1.5); end
        if plotInferred  plot(t,inferred(idx,:)-idx*sep,'color',[0 0 1],'linewidth',1); end
    end    
    ylabel('\DeltaF/F');
end
xlabel('Time (s)');

if displayLabel
    for idx=1:K
        text(max(t)*1.02,-idx*sep, num2str(idx));
    end
end

h=gca;
if ~rollingView 
    axis(h,[0 max(t) -sep*(K+1) sep]);
else
    idx=1;
    displayNum=20;
    while idx+displayNum<K
        axis(h,[0 max(t) -sep*(idx+displayNum) sep*(-idx+1)]);
        idx=idx+displayNum;
        input('Press ENTER to continue...');
    end
    axis(h,[0 max(t) -sep*(idx+displayNum) sep*(-idx+1)]);
    input('Press ENTER to continue...');
    axis(h,[0 max(t) -sep*(K+1) sep]);    
end

