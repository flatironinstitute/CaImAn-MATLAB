function figH = plot4Dproj(Y, maskY, sizY, params)
% plot4Dproj(Y, maskY, sizY, params)
% displays projections of 3D-4D matrices (Y), if 4D the figure will update with
% a lag of pi.lag and overlais a maskY
global pi
pi = [];
pi.iter = 0;
pi.range = [];
pi.cbgate = 1;
pi.lag = 0.01;
pi.maskcor = [1 0 0];
pi.cormap = 'jet';
pi.maskout = 0;

fop = fields(pi);
if exist('params', 'var') && ~isempty(params)
    for ii = 1:length(fop)
        if isfield(params, fop{ii}) && ~isempty(params(1). (fop{ii}))
            pi.(fop{ii}) = params(1).(fop{ii});
        end
    end
end

% flatten all
Y = reshape(full(Y), prod(sizY), []);
maskY = reshape(full(maskY), prod(sizY), []);

% reshape Y and maskY if it is flatten
Y = reshape(full(Y), [sizY, size(Y, 2)]);
maskY = reshape(full(maskY), [sizY, size(maskY, 2)]);

Y = double(Y);

%% Plotting
figH = figure('position', [368 286 1156 510]); % [213 311 1461 787]);
colormap(pi.cormap)

if pi.maskout && size(maskY, 4) == 1
   Y = bsxfun(@times, Y, double(maskY > 0));
end

if isempty(pi.range) && pi.iter == 0
   pi.range = [min(Y(:)) max(Y(:))];
end

plotmuliproj(Y, maskY)
end

function plotmuliproj(Y, maskY)
global pi
pi.d2proj = [3 2 1];
pi.axesname = {'X', 'Y'; 'Z', 'Y'; 'Z', 'X'};
pi.c2plot = [2, 1; 3, 1; 3, 2];
if isempty(pi.range)
    pi.range = [min(Y(:)) max(Y(:))];
end
hAxes(1) = subplot(1, 3, 1);
hAxes(2) = subplot(1, 3, 2);
hAxes(3) = subplot(1, 3, 3);
fprintf(['Plotting ', num2str(size(Y, 4)), ' volume(s)\n'])
for t = 1:size(Y, 4)
    for a_idx = 1:3
        %% plot background Y
        plotproj(hAxes, Y(:, :, :, t), a_idx, t)
        %% plot overlay maskY
        if ~isempty(maskY)
            if size(Y, 4) == size(maskY, 4) % if mask is same size uses each for each Y(:, :, :, t)
                plotoverlay(hAxes, a_idx, maskY(:, :, :, t))
            elseif size(maskY, 4) == 1 % otherwise it uses the same for all (1)
                plotoverlay(hAxes, a_idx, maskY(:, :, :))
            end
        end
    end
    pause(pi.lag)
end
end

function plotproj(hAxes, Im, hi, ti)
global pi
imagesc(squeeze(max(Im, [], pi.d2proj(hi))), 'Parent', hAxes(hi))
if pi.iter; pi.range = [min(Im(:)) max(Im(:))]; end % update caxis iteratively at each t
if ~isempty(pi.range); caxis(hAxes(hi), pi.range); end
if pi.cbgate; colorbar(hAxes(hi)); end
xlabel(hAxes(hi), pi.axesname{hi, 1}); 
ylabel(hAxes(hi), pi.axesname{hi, 2});
title(hAxes(hi), num2str(ti));
set(hAxes(hi), 'YTick', []); set(hAxes(hi), 'XTick', [])
end

function plotoverlay(hAxes, idx, maskY)
global pi
hold(hAxes(idx) , 'on')
contour(squeeze(max(maskY, [], pi.d2proj(idx))) > 0, ...
   1, 'color', pi.maskcor, 'Parent', hAxes(idx))
end