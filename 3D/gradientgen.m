function colormatrix = gradientgen(cmap2use, ncolors)
if ~exist('cmap2use', 'var') || isempty(cmap2use); cmap2use = []; end
if ~exist('ncolors', 'var') || isempty(ncolors); ncolors = []; end
% uses buildcmap to generate a gradient between given colors
% buildcmap input can be str or vector values
cmaps = {'wr', 'wk', 'wb', 'wg', 'wm', 'wc', 'kr', 'mbg', 'gry'}; % but you can even mix 3 of them
maxc = 180;
if ~isempty(cmap2use)
    if isstr(cmap2use)
        cmap2use = cmap2use;
    else
        cmap2use = cmaps{cmap2use};
    end
    cmap = buildcmap(cmap2use); 
    cmap = cmap(round(76:255), :);
    if cmap2use == 1
        colormatrix = cmap(maxc,:);
    else % resampling
        colormatrix = interp1(1:maxc,cmap,1:maxc/(ncolors+1):maxc);
        colormatrix = colormatrix(1:ncolors,:);
    end
else
    fprintf('Posible colormaps (')
    for i = 1:numel(cmaps); if i ~= numel(cmaps); fprintf([cmaps{i}, ' , ']); else fprintf([cmaps{i}]); end; end
    fprintf(') \n')
    colormatrix = [];
end