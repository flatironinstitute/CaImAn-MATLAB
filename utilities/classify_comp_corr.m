function [rval_space,rval_time,ind_space,ind_time] = classify_comp_corr(Yr,A,C,b,f,options)

% Classify components by looking at how well the spatial and temporal components
% correlate with the raw data at point of high activation

% INPUTS
% Yr:       reshaped data in 2D matrix or memmaped file
% A,C,b,f:  identified components
% options:  options structure
%    space_thresh:  r-value threshold for spatial components (default: 0.4) 
%    time_thresh:  r-value threshold for temporal components (default: 0.4) 
%    Athresh:  threshold for determining spatial overlap (default: 0.1)
%    Np:       number of high activity intervals for each component (default: 20)
%    peak_int: Interval around each local peak to be considered (default: -2:6)
%    MinPeakDist:   minimum peak distance for finding points of high activity  (default: 10)

% OUTPUTS:
% rval_space:     r-values between spatial components and data patches
% rval_time:      r-values between temporal components and trace averages
% ind_space:      components with rval_space > space_thresh
% ind_time:       components with rval_time > time_thresh
% MinPeakDist:

% Written by Eftychios A. Pnevmatikakis, Simons Foundation, 2016
% based on discussions with Matt Kaufman and Farzaneh Najafi, CSHL

defoptions = CNMFSetParms;
if nargin < 6 || isempty(options)
    options = defoptions;
end

if ~isfield(options,'peak_int') || isempty(options.peak_int); int = defoptions.peak_int; else int = options.peak_int; end
if ~isfield(options,'Npeaks') || isempty(options.peak_int); Np = defoptions.Npeaks; else Np = options.Npeaks; end
if ~isfield(options,'A_thresh') || isempty(options.A_thresh); Athresh = defoptions.A_thresh; else Athresh = options.A_thresh; end
if ~isfield(options,'space_thresh') || isempty(options.space_thresh); options.space_thresh = defoptions.space_thresh; end
if ~isfield(options,'time_thresh') || isempty(options.time_thresh); options.time_thresh = defoptions.time_thresh; end
if ~isfield(options,'MinPeakDist') || isempty(options.MinPeakDist); options.MinPeakDist = defoptions.MinPeakDist; end

memmaped = isobject(Yr);

[K_m,T] = size(C);
pk = zeros(K_m,Np);

parfor i = 1:K_m
    [~,pk_temp] = findpeaks(C(i,:),'SortStr','descend','Npeaks',Np,'MinPeakDistance',options.MinPeakDist);
    if ~isempty(pk_temp)
        if length(pk_temp) < Np
            pk_temp(length(pk_temp)+1:Np) = pk_temp(1);
        end
        pk(i,:) = pk_temp;
    end
end

%% expand intervals
lint = length(int);

locs = repmat(pk,1,lint) + kron(int,ones(size(pk)));
locs = sort(locs,2,'ascend');
locs(locs<1)=1;
locs(locs>T)=T;
LOCS = mat2cell(locs,ones(K_m,1),Np*lint);
for i = 1:K_m
    LOCS{i} = unique(LOCS{i});
end

%% compute r-values
nA = full(sqrt(sum(A.^2)))';
tAA = (A'*A)./(nA*nA') > Athresh;
tAA(1:K_m+1:K_m^2) = 0;

rval_space = zeros(K_m,1);
rval_time  = zeros(K_m,1);
if memmaped
    data = Yr.Yr;
else
    data = Yr;
end
for i = 1:K_m
    ovlp_cmp = find(tAA(:,i));
    indeces = LOCS{i};
    for j = 1:length(ovlp_cmp)
        indeces = setdiff(indeces,LOCS{ovlp_cmp(j)});
    end
    amask = A(:,i)>0;
    if ~isempty(indeces)
        if memmaped      
    %         time_indeces = sort(indeces,'ascend');
    %         ff_time = [0,find(diff(time_indeces)>1),length(time_indeces)];
    %         time_intervals = mat2cell(time_indeces,1,diff(ff_time));
    %         space_indeces = sort(find(amask),'ascend')';
    %         ff_space = [0,find(diff(space_indeces)>1),length(space_indeces)];
    %         space_intervals = mat2cell(space_indeces,1,diff(ff_space));
    %         data = zeros(ff_space(end),ff_time(end));
    %         for rows = 1:length(space_intervals)
    %             for columns = 1:length(time_intervals)
    %                 data(ff_space(rows)+1:ff_space(rows+1),ff_time(columns)+1:ff_time(columns+1)) = Yr.Yr(space_intervals{rows},time_intervals{columns});
    %             end
    %         end                
            %mY = (A(:,i)>0).*double(mean(Yr.Yr(:,indeces),2) - b*mean(f(:,indeces),2));
            mY_space = double(mean(data(amask,indeces),2) - b(amask,:)*mean(f(:,indeces),2));
            mY_time = double(mean(data(amask,indeces)) - mean(b(amask,:))*f(:,indeces));
        else
            mY_space = double(mean(Yr(amask,indeces),2) - b(amask,:)*mean(f(:,indeces),2));
            mY_time = double(mean(Yr(amask,indeces)) - mean(b(amask,:))*f(:,indeces));
        end
        rval_space(i) = corr(full(A(amask,i)),mY_space);
        rval_time(i) = corr(C(i,indeces)',mY_time');
    else
        rval_space(i) = NaN;
        rval_time(i) = NaN;
    end
end

ind_space = rval_space > options.space_thresh;
ind_time = rval_time > options.time_thresh;
