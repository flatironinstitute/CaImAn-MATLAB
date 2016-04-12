function [C_df,Df,S_df] = extract_DF_F(Y,A,C,S,i,options)

% extract DF/F signals after performing NMF
% inputs:  Y raw data (d X T matrix, d # number of pixels, T # of timesteps)
%          A matrix of spatial components (d x K matrix, K # of components)
%          C matrix of temporal components (K x T matrix)
%          S matrix of deconvolved activity ((K-1) x T matrix) (optional)
%          i index of component that represent the background (optional, if not
%          given it's estimated)
%          options structure used for specifying method for determining DF
%           default method is the median of the trace. By changing
%           options.df_prctile an arbitray percentile can be used (between 0 and 100).
%           a moving window can also be established by specifying options.df_window

% outputs:  C_df temporal components in the DF/F domain
%           Df   background for each component to normalize the filtered raw data    
%           S_df deconvolved activity/spikes in the DF/F domain

% Written by: 
% Eftychios A. Pnevmatikakis, Simons Foundation, 2015

defoptions = CNMFSetParms;
if nargin < 6 || isempty(options)
    options = defoptions;
end
if ~isfield(options,'df_prctile') || isempty(options.df_prctile)
    options.df_prctile = defoptions.df_prctile;
end
if ~isfield(options,'df_window') || isempty(options.df_window)
    options.df_window = defoptions.df_window;
end

nA = sqrt(sum(A.^2))';
[K,~] = size(C);
A = A/spdiags(nA,0,K,K);    % normalize spatial components to unit energy
C = spdiags(nA,0,K,K)*C;

if nargin < 5 || isempty(i)
    [~,i] = min(sum(A.^6)); % identify background component
end

non_bg = true(1,K); 
non_bg(i) = false;      % non-background components
Yf = A'*Y - (A'*A(:,non_bg))*C(non_bg,:);

if isempty(options.df_window) || (options.df_window > size(C,2))
    if options.df_prctile == 50
        Df = median(Yf,2);
    else
        Df = prctile(Yf,options.df_prctile,2);
    end
    C_df = spdiags(Df,0,K,K)\C;
else
    if options.df_prctile == 50
        Df = medfilt1(Yf,options.df_window,[],2,'truncate');
    else
        Df = zeros(size(Yf));
        for i = 1:size(Df,1);
            df_temp = running_percentile(Yf(i,:), options.df_window, options.df_prctile);
            Df(i,:) = df_temp(:)';
        end
    end
    C_df = C./Df;
end
            
C_df(i,:) = 0;

if nargin < 4 || isempty(S) || nargout < 3
    S_df = [];
    if nargout == 3
        warning('Merged spikes matrix is returned as empty because the original matrix was not provided.');
    end
else
    if isempty(options.df_window) || (options.df_window > size(C,2))
        S_df = spdiags(Df(non_bg(:)),0,sum(non_bg),sum(non_bg))\S;
    else
        S_df = S./Df(non_bg,:);
    end
end