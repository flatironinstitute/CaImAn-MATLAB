function [A,C,nr,merged_ROIs,P,S] = merge_components(Y,A,b,C,f,P,S,options)

% merging of spatially overlapping components that have highly correlated tmeporal activity
% The correlation threshold for merging overlapping components is user specified in P.merge_thr (default value 0.85)
% Inputs:
% Y:            raw data
% A:            matrix of spatial components
% b:            spatial background
% C:            matrix of temporal components
% f:            temporal background
% P:            struct for neuron parameters
% S:            deconvolved activity/spikes (optional)
% options:      struct for algorithm parameters

% Outputs:
% A:            matrix of new spatial components
% C:            matrix of new temporal components
% nr:           new number of components
% merged_ROIs:  list of old components that were merged
% P:            new parameter struct
% S:            matrix of new deconvolved/activity spikes

% Written by:
% Eftychios A. Pnevmatikakis, Simons Foundation, 2015

defoptions = CNMFSetParms;
if nargin < 8; options = []; end
if ~isfield(options,'d1') || isempty(options.d1); d1 = input('What is the total number of rows? \n'); else d1 = options.d1; end          % # of rows
if ~isfield(options,'d2') || isempty(options.d2); d2 = input('What is the total number of columns? \n'); else d2 = options.d2; end       % # of columns
if ~isfield(options,'merge_thr') || isempty(options.merge_thr); thr = defoptions.merge_thr; else thr = options.merge_thr; end     % merging threshold
if ~isfield(options,'max_merg'); mx = 50; else mx = options.max_merg; end           % maximum merging operations
if ~isfield(options,'deconv_method') || isempty(options.deconv_method); options.deconv_method = defoptions.deconv_method; end
if ~isfield(options,'fast_merge') || isempty(options.fast_merge); options.fast_merge = defoptions.fast_merge; end  % flag for using fast merging

nr = size(A,2);
%[d,T] = size(Y);
d = size(A,1);
T = size(C,2);
C_corr = corr(full(C(1:nr,:)'));
FF1 = triu(C_corr)>= thr;                           % find graph of strongly correlated temporal components

A_corr = triu(A(:,1:nr)'*A(:,1:nr));                
A_corr(1:nr+1:nr^2) = 0;
FF2 = A_corr > 0;                                   % find graph of overlapping spatial components

FF3 = and(FF1,FF2);                                 % intersect the two graphs
[l,c] = graph_connected_comp(sparse(FF3+FF3'));     % extract connected components
MC = [];
for i = 1:c
    if length(find(l==i))>1
        MC = [MC,(l==i)'];
    end
end

cor = zeros(size(MC,2),1);
for i = 1:length(cor)
    fm = find(MC(:,i));
    for j1 = 1:length(fm)
        for j2 = j1+1:length(fm)
            cor(i) = cor(i) + C_corr(fm(j1),fm(j2));
        end
    end
end

[~,ind] = sort(cor,'descend');
nm = min(length(ind),mx);   % number of merging operations
merged_ROIs = cell(nm,1);
A_merged = zeros(d,nm);
C_merged = zeros(nm,T);
S_merged = zeros(nm,T);
if strcmpi(options.deconv_method,'constrained_foopsi')
    P_merged.gn = cell(nm,1);
    P_merged.b = cell(nm,1);
    P_merged.c1 = cell(nm,1);
    P_merged.neuron_sn = cell(nm,1);
end
if ~options.fast_merge
    Y_res = Y - A*C;
end

for i = 1:nm
    merged_ROIs{i} = find(MC(:,ind(i)));
    nC = sqrt(sum(C(merged_ROIs{i},:).^2,2));
    if options.fast_merge
        aa = sum(A(:,merged_ROIs{i})*spdiags(nC,0,length(nC),length(nC)),2);
        for iter = 1:10
            cc = (aa'*A(:,merged_ROIs{i}))*C(merged_ROIs{i},:)/sum(aa.^2);
            aa = A(:,merged_ROIs{i})*(C(merged_ROIs{i},:)*cc')/norm(cc)^2;
        end
        na = sqrt(sum(aa.^2));
        aa = aa/na;
        %[cc,b_temp,c1_temp,g_temp,sn_temp,ss] = constrained_foopsi(cc);
        cc = na*cc';
        ss = cc;
    else
        A_merged(:,i) = sum(A(:,merged_ROIs{i})*spdiags(nC,0,length(nC),length(nC)),2);    
        Y_res = Y_res + A(:,merged_ROIs{i})*C(merged_ROIs{i},:);
        ff = find(A_merged(:,i));
        Pmr = P;
        if isfield(Pmr,'unsaturatedPix');
            px = intersect(Pmr.unsaturatedPix,ff);
            Pmr.unsaturatedPix = zeros(length(px),1);
            for pxi = 1:length(px)
                Pmr.unsaturatedPix(pxi) = find(ff == px(pxi));
            end
        end
        cc = update_temporal_components(Y_res(ff,:),A_merged(ff,i),b(ff,:),median(spdiags(nC,0,length(nC),length(nC))\C(merged_ROIs{i},:)),f,Pmr,options);
        [aa,bb] = update_spatial_components(Y_res,cc,f,A_merged(:,i),P,options);    
        [cc,~,Ptemp,ss] = update_temporal_components(Y_res(ff,:),aa(ff),bb(ff,:),cc,f,Pmr,options);
    end
    A_merged(:,i) = aa;    
    C_merged(i,:) = cc;
    S_merged(i,:) = ss;
    if strcmpi(options.deconv_method,'constrained_foopsi') || strcmpi(options.deconv_method,'MCEM_foopsi')
        if options.fast_merge
            P_merged.gn{i} = 0; %g_temp;   % do not perform deconvolution during merging
            P_merged.b{i} = 0;  %b_temp;
            P_merged.c1{i} = 0; %c1_temp;
            P_merged.neuron_sn{i} = 0; %sn_temp;
        else
            P_merged.gn{i} = Ptemp.gn{1};
            P_merged.b{i} = Ptemp.b{1};
            P_merged.c1{i} = Ptemp.c1{1};
            P_merged.neuron_sn{i} = Ptemp.neuron_sn{1};
            if i < nm
                Y_res(ff,:) = Y_res(ff,:) - aa(ff)*cc;
            end
        end
    end
end

neur_id = unique(cell2mat(merged_ROIs));

A = [A(:,1:nr),A_merged,A(:,nr+1:end)];
C = [C(1:nr,:);C_merged;C(nr+1:end,:)];
A(:,neur_id) = [];
C(neur_id,:) = [];

if nargin < 7
    S = [];
    if nargout == 6
        warning('Merged spikes matrix is returned as empty because the original matrix was not provided.');
    end
else
    S = [S(1:nr,:);S_merged];
    S(neur_id,:) = [];
end

if strcmpi(options.deconv_method,'constrained_foopsi') || strcmpi(options.deconv_method,'MCEM_foopsi')
    P.b(neur_id) = [];
    P.b(nr - length(neur_id) + (1:nm)) = P_merged.b;
    P.gn(neur_id) = [];
    P.gn(nr - length(neur_id) + (1:nm)) = P_merged.gn;
    P.c1(neur_id) = [];
    P.c1(nr - length(neur_id) + (1:nm)) = P_merged.c1;
    P.neuron_sn(neur_id) = [];
    P.neuron_sn(nr - length(neur_id) + (1:nm)) = P_merged.neuron_sn;
end
nr = nr - length(neur_id) + nm;