function figH = plotCenteroverY(Y, icenter, sizY, params)
% plotCenteroverY(Y, icenter, sizY, params)
% Displays center or roi's over a XY, YZ and XZ projections of Y.
% inputs:
% Y : background image 2D or 3D matrix (spatial dimensions are 2D or 3D)
% icenter: center or rois kxm (k == number of rois, m == spatial dimensions (2D or 3D))
% sizY: spatial dimensions of Y or icenter (2D or 3D);
% params: extra parameters that edit the figure fitures

global pi
pi = [];
pi.range = []; % fluorescence range (for caxis)
pi.ititle = 'ROI centers'; % time lag
pi.d2proj = [3 2 1]; % dimension to project
pi.axesname = {'X', 'Y'; 'Z', 'Y'; 'Z', 'X'}; % axes name
if length(sizY) == 3
    pi.rColor = gradientgen(9, sizY(3)); % color of rois: colorcode the depth
else
    pi.rColor = gradientgen(9, size(icenter, 1)); % color of rois: colorcode the roi idx
end
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

for a_idx = 1:length(hAxes)
    plotproj(hAxes, Y, a_idx)
    plotcenter(hAxes, icenter, a_idx)
end
end

function plotproj(hAxes, Y, idx)
global pi
%% Plotting
imagesc(squeeze(max(Y, [], pi.d2proj(idx))), 'Parent', hAxes(idx))
caxis(hAxes(idx), pi.range)
%% Figure details
if ~isempty(pi.range); caxis(hAxes(idx), pi.range); end
xlabel(hAxes(idx), pi.axesname{idx, 1}); 
ylabel(hAxes(idx), pi.axesname{idx, 2});
title(hAxes(idx), pi.ititle);
set(hAxes(idx), 'YTick', []); set(hAxes(idx), 'XTick', [])
end

function plotcenter(hAxes, icenter, idx)
global pi
hold(hAxes(idx), 'on')
scatter(icenter(:, pi.c2plot(idx, 1)), icenter(:, pi.c2plot(idx, 2)), ...
50, pi.rColor(round(icenter(:, 3)), :), 'Parent', hAxes(idx));
end