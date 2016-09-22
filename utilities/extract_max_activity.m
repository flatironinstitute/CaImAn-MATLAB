function [LOCS,Ym] = extract_max_activity(Y,nlocs,width)

% extract points of maximum activation for each pixel using the findpeaks
% function. Requires the signal processing toolbox and uses the parallel
% processing toolbox if present

% INPUTS:
% Y:        input data
% nlocs:    number of local maxima for each location
% width:    length of each interval

% OUTPUTS:
% LOCS:     max activation intervals
% Ym:       max activation values

% Written by:
% Eftychios A. Pnevmatikakis, Simons Foundation, 2016

warning('off', 'MATLAB:maxNumCompThreads:Deprecated');

if nargin < 3 || isempty(width)
    width = 21;
end

if nargin < 2 || isempty(nlocs)
    nlocs = 30;
end

sizY = size(Y);
T = sizY(end);
d = prod(sizY(1:end-1));
if length(sizY) > 2
    Y = reshape(Y,d,T);
end

if ~mod(width,2); width = width + 1; end
locs = zeros(d,nlocs);

spatial_parallel = exist(which('parpool'));

if spatial_parallel         % solve BPDN problem for each pixel
    Nthr = max(20*maxNumCompThreads,round(d*T/2^24));
    Nthr = min(Nthr,round(d/1e3));
    siz_row = [floor(d/Nthr)*ones(Nthr-mod(d,Nthr),1);(floor(d/Nthr)+1)*ones(mod(d,Nthr),1)];
    indeces = [0;cumsum(siz_row)];    
    loc_cell = cell(Nthr,1);
    for nthr = 1:Nthr 
        lc_temp = zeros(siz_row(nthr),nlocs);
        Ytemp = double(Y(indeces(nthr)+1:indeces(nthr+1),:));
        parfor px = 1:siz_row(nthr)
            [pks,lcs] = findpeaks(Ytemp(px,:),'MinPeakDistance',width);
            [~,srt] = sort(pks,'descend');
            if length(srt) >= nlocs
                lc_temp(px,:) = lcs(srt(1:nlocs));
            else
                rnd = randperm(floor(T/width)-1,nlocs-length(srt))*width;
                lc_temp(px,:) = [srt,rnd];
            end            
        end
        loc_cell{nthr} = lc_temp;
        if mod(nthr,50) == 0
            fprintf('%2.1f%% of pixels completed \n', indeces(nthr+1)*100/d);
        end
    end
    locs = cell2mat(loc_cell);
else
    Y = double(Y);
    for i = 1:d
        [pks,lcs] = findpeaks(Y(i,:),'MinPeakDistance',width);
        [~,srt] = sort(pks,'descend');
        if length(srt) >= nlocs
            locs(i,:) = lcs(srt(1:nlocs));
        else
            rnd = randperm(floor(T/width)-1,nlocs-length(srt))*width;
            locs(i,:) = [srt,rnd];
        end
    end
end

locs(locs<(width+1)/2) = (width+1)/2;
locs(locs>T - (width+1)/2) = T - (width+1)/2;
locs = sort(locs,2,'ascend');
LOCS = kron(locs,ones(1,width)) + repmat(kron(ones(1,nlocs),-(width-1)/2:(width-1)/2),d,1);
if nargout == 2
    Ym = zeros(size(LOCS));
    for i = 1:d
        Ym(i,:) = Y(i,LOCS(i,:));
    end
end