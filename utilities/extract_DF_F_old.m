function [C_df,Df] = extract_DF_F_old(Y,A,C,ind,options)

% extract DF/F signals after performing NMF
% inputs:  Y raw data (d X T matrix, d # number of pixels, T # of timesteps)
%          A matrix of spatial components (d x K matrix, K # of components)
%          C matrix of temporal components (K x T matrix)
%          ind index of component that represent the background (optional, if not
%          given it's estimated)
%          options structure used for specifying method for determining DF
%           default method is the median of the trace. By changing
%           options.df_prctile an arbitray percentile can be used (between 0 and 100).
%           a moving window can also be established by specifying options.df_window

% outputs:  C_df temporal components in the DF/F domain
%           Df   background for each component to normalize the filtered raw data    

% Written by: 
% Eftychios A. Pnevmatikakis, Simons Foundation, 2015

memmaped = isobject(Y);
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
if ~isfield(options,'full_A') || isempty(options.full_A); full_A = defoptions.full_A; else full_A = options.full_A; end

nA = sqrt(sum(A.^2))';
[K,T] = size(C);
d = size(A,1);
A = A/spdiags(nA,0,K,K);    % normalize spatial components to unit energy
C = spdiags(nA,0,K,K)*C;

if nargin < 4 || isempty(ind)
    [~,ind] = min(sum(A.^6)); % identify background component
end

non_bg = true(1,K); 
non_bg(ind) = false;      % non-background components

step = 5e3;
if memmaped
    AY = zeros(K,T);
    for i = 1:step:d
        AY = AY + A(i:min(i+step-1,d),:)'*double(Y.Yr(i:min(i+step-1,d),:));
    end
else
    if issparse(A) && isa(Y,'single')  
        if full_A
            AY = full(A)'*Y;
        else
            AY = A'*double(Y);
        end
    else
        AY = A'*Y;
    end
end

Yf = AY - (A'*A(:,non_bg))*C(non_bg,:);

if isempty(options.df_window) || (options.df_window > size(C,2))
    if options.df_prctile == 50
        Df = median(Yf,2);
    else
        Df = prctile(Yf,options.df_prctile,2);
    end
    C_df = spdiags(Df,0,K,K)\C;
else
    if options.df_prctile == 50
        if verLessThan('matlab','2015b')
            warning('Median filtering at the boundaries might be inaccurate due to zero padding.')
            Df = medfilt1(Yf,options.df_window,[],2);
        else
            Df = medfilt1(Yf,options.df_window,[],2,'truncate');
        end
    else
        Df = zeros(size(Yf));
        for i = 1:size(Df,1);
            df_temp = running_percentile(Yf(i,:), options.df_window, options.df_prctile);
            Df(i,:) = df_temp(:)';
        end
    end
    C_df = C./Df;
end
            
C_df(ind,:) = []; % 0; % FN modified so C_df does not include the background components and it has the same size as C.