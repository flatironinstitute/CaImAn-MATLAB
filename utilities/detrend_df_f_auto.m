function [F_dff,F0] = detrend_df_f_auto(A,b,C,f,YrA,options)

% detrend and extract DF/F values from processed fluorescence traces
% INPUTS
% A:        matrix of spatial components (2-d matrix, sparse)
% b:        matrix of spatial background (2-d matrix)
% C:        matrix of temporal components (2-d matrix)
% f:        matrix of temporal background (2-d matrix)
% YrA:      matrix of filtered residuals (2-d matrix, optional)
% options:  options structure used for specifying method for determining DF
%           default method is the median of the trace. By changing
%           options.df_prctile an arbitray percentile can be used (between 0 and 100).
%           a moving window can also be established by specifying options.df_window
%
% OUTPUTS
% F_dff:    detrended fluorescence in DF/F 
% F0:       baseline fluorescence for each trace

defoptions = CNMFSetParms;
if nargin < 5 || isempty(options)
    options = defoptions;
end
if ~isfield(options,'df_prctile') || isempty(options.df_prctile)
    options.df_prctile = defoptions.df_prctile;
end
if ~isfield(options,'df_window') || isempty(options.df_window)
    options.df_window = defoptions.df_window;
end

F = diag(sum(A.^2))*(C + YrA);                                                     % fluorescence
B = (A'*b)*f;

[N,T] = size(F);

F_dff = zeros(N,T);
F0 = zeros(N,T);

for i = 1:N
    if isempty(options.df_window) || (options.df_window > size(C,2))
        [~,cdf] = estimate_percentile_level(F(i,:),T,T);
        cdf_level = median(cdf);
        Fd = prctile(F(i,:),cdf_level,2);
        F0(i,:) = repmat(prctile(B(i,:),cdf_level,2) + Fd,1,T);
        F_dff(i,:) = (F(i,:) - Fd)./F0(i,:);
    else
        [~,cdf] = estimate_percentile_level(F(i,:),options.df_window,round(options.df_window/2));
        cdf_level = median(cdf);
        Fd = prctfilt(F(i,:),cdf_level,options.df_window);                               % detrended fluorescence
        F0(i,:) = prctfilt(B(i,:),cdf_level,options.df_window,[],0) + (F(i,:)-Fd);       % background + baseline for each component
        F_dff(i,:) = Fd./F0(i,:);
    end
end