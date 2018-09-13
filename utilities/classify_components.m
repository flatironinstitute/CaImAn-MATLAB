function [rval_space,rval_time,max_pr,sizeA,keep] = classify_components(Y,A,C,b,f,YrA,options)

% perform several tests for component classification:
%  i)   correlation test using the function classify_comp_corr
%  ii)  calculate the probability that max of each temporal component is
%         due to noise through trace_fit_extreme.m
%  iii) filter component based on min/max size

% OUTPUTS
% rval_space:       correlation in space
% rval_time:        correlation in time
% max_pr:           max probability
% sizeA:            size of each component
% keep:             binary flag vector on whether each component passes all tests.


defoptions = CNMFSetParms;
if nargin < 7 || isempty(options)
    options = defoptions;
end

if ~isfield(options,'max_pr_thr') || isempty(options.max_pr_thr); options.max_pr_thr = defoptions.max_pr_thr; end
if ~isfield(options,'fr') || isempty(options.fr); fr = defoptions.fr; else; fr = options.fr; end
if ~isfield(options,'t_int') || isempty(options.t_int); t_int = defoptions.t_int; else; t_int = options.t_int; end
if ~isfield(options,'sn_fac') || isempty(options.t_int); sn_fac = defoptions.sn_fac; else; sn_fac = options.sn_fac; end
if ~isfield(options,'max_size_thr') || isempty(options.max_size_thr); options.max_size_thr = defoptions.max_size_thr; end
if ~isfield(options,'min_size_thr') || isempty(options.min_size_thr); options.min_size_thr = defoptions.min_size_thr; end
if ~isfield(options,'size_thr') || isempty(options.size_thr); options.size_thr = defoptions.size_thr; end

[rval_space,rval_time,ind_space,ind_time] = classify_comp_corr(Y,A,C,b,f,options);

%AA = A'*A;
%AY = mm_fun(A,Y);
% if nargin < 6 || isempty(YrA)
%     YrA = bsxfun(@times, 1./sum(A.^2)',AY - AA*C);
% end

%traces = detrend_df_f(A,b,C,f,YrA,options);
%max_pr = trace_fit_extreme(full(traces),fr,t_int,sn_fac);
max_pr = ones(size(rval_space));
%fitness = compute_event_exceptionality(traces,0);
%fitness_delta = compute_event_exceptionality(diff(traces,[],2),0);

sizeA = zeros(size(A,2),1);
for i = 1:length(sizeA)
    sizeA(i) = sum(A(:,i)>options.size_thr*max(A(:,i)));
end

keep = ind_space & ind_time & (max_pr > options.max_pr_thr) & (sizeA >= options.min_size_thr) & (sizeA <= options.max_size_thr);