function plotCenteroverY(Y, icenter, sizY, params)
% plotCenteroverY(Y, icenter, sizY, params)
% This function plots the center of given ROIs (2D or 3D) icenter over the max
% projections of a background image Y
global pi
pi = [];
pi.contcolor = {'y', 'g', 'r'}; % contour color
pi.range = []; % fluorescence range (for caxis)
pi.lag = 0.2; % time lag
pi.d2proj = [3 2 1]; % dimension to project
pi.axesname = {'X', 'Y'; 'Z', 'Y'; 'Z', 'X'}; % axes name
pi.rColor = gradientgen(9, sizY(3)); % color of rois: colorcode the depth
pi.c2plot = [2, 1; 3, 1; 3, 2];
pi.figpos = [-1267 326 1156 510]; % figure position

fop = fields(pi);
if exist('params', 'var') && ~isempty(params)
    for ii = 1:length(fop)
        if isfield(params, fop{ii}) && ~isempty(params(1). (fop{ii}))
            pi.(fop{ii}) = params(1).(fop{ii});
        end
    end
end

%% Plotting
figH = figure('position', pi.figpos);
colormap('gray')

%% always reshape Y to target sizY
Y = double(reshape(Y, sizY));

plotmuliproj(Y, icenter, sizY)
end

function plotmuliproj(Y, icenter, sizY)
global pi
%% getting number of components
tn = size(icenter, numel(sizY) + 1);

%% figure features
if numel(sizY) == 3
    hAxes(1) = subplot(1, 3, 1);
    hAxes(2) = subplot(1, 3, 2);
    hAxes(3) = subplot(1, 3, 3);
else
    hAxes(1) = subplot(1, 2, 1);
end

if isempty(pi.range)
    pi.range = [min(Y(:)) max(Y(:))];
    display(pi.range)
end

for t = 1:tn
    for a_idx = 1:length(hAxes)
        plotproj(hAxes, Y, a_idx, num2str(t))
        plotcenter(hAxes, icenter, a_idx)
    end
    pause(pi.lag)
end
end

function plotproj(hAxes, Y, idx, tidx)
global pi
%% Plotting
imagesc(squeeze(max(Y, [], pi.d2proj(idx))), 'Parent', hAxes(idx))
caxis(hAxes(idx), pi.range)
%% Figure details
if ~isempty(pi.range); caxis(hAxes(idx), pi.range); end
xlabel(hAxes(idx), pi.axesname{idx, 1}); 
ylabel(hAxes(idx), pi.axesname{idx, 2});
title(hAxes(idx), tidx);
set(hAxes(idx), 'YTick', []); set(hAxes(idx), 'XTick', [])
end

function plotcenter(hAxes, icenter, idx)
global pi
hold(hAxes(idx), 'on')
scatter(icenter(:, pi.c2plot(idx, 1)), icenter(:, pi.c2plot(idx, 2)), ...
50, pi.rColor(round(icenter(:, 3)), :), 'Parent', hAxes(idx));
end