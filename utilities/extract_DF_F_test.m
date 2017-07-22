function [C_df,Df] = extract_DF_F_test(Y,A,C,b,f,YrA,P,options)

% extract DF/F signals after performing NMF
% inputs:  Y:    raw data (d X T matrix, d # number of pixels, T # of timesteps)
%          A:    matrix of spatial components (d x K matrix, K # of components)
%          C:    matrix of temporal components (K x T matrix)
%          b:    spatial background
%          YrA:  residual fluorescence for each component   
%          P:   neuron structure, used to read the baseline activity for each component of C
%          options structure used for specifying method for determining DF
%               default method is the median of the trace. By changing
%               options.df_prctile an arbitray percentile can be used (between 0 and 100).
%               a moving window can also be established by specifying options.df_window
%               options.add_residual specifies whether to add the residual
%               fluoresence YrA to the baseline computation or not (default: true)

% outputs:  C_df temporal components in the DF/F domain
%           Df   background for each component to normalize the filtered raw data    

% Written by: 
% Eftychios A. Pnevmatikakis, Simons Foundation, 2016

%memmaped = isobject(Y);
defoptions = CNMFSetParms;
if ~exist('options','var') || isempty(options); options = defoptions; end
if ~isfield(options,'df_prctile') || isempty(options.df_prctile); options.df_prctile = defoptions.df_prctile; end
if ~isfield(options,'df_window') || isempty(options.df_window); options.df_window = defoptions.df_window; end
if ~isfield(options,'add_residual') || isempty(options.add_residual); options.add_residual = defoptions.add_residual; end

[K,T] = size(C);
options.df_window = min(options.df_window,T);
nA = sum(A.^2);
AA = A'*A;

if (~exist('YrA','var') || isempty(YrA)) && options.add_residual    % compute and add residual fluorescence
    AY = mm_fun(A,Y);
    YrA = spdiags(nA(:),0,K,K)\(AY - AA*C - (A'*b)*f);
else
    YrA = sparse(K,T);
end

Bas = zeros(K,T);
bas_val = zeros(K,1);

if ~(~exist('P','var') || isempty(P))
    bas_val = cell2mat(P.b);
    Ntr = size(bas_val,2);
    if Ntr > 1
        ln = diff(P.cs_frtrs);
        for i = 1:Ntr
            Bas(:,:,P.cs_frtrs(i)+1:P.cs_frtrs(i+1)) = repmat(bas_val(:,i),1,ln(i));
        end
    else
        Bas = repmat(bas_val,1,T);
    end
end

Df = C - Bas + options.add_residual*YrA;    % DF signal
%AAov = AA - diag(diag(AA)); % matrix of overlapping components
C_bg = bsxfun(@times,1./nA(:),(A'*b)*f) + Bas + options.add_residual*YrA;  % background fluorescence for each component

if isempty(options.df_window) || (options.df_window >= T)
    if options.df_prctile == 50
        F0 = median(C_bg,2);
    else
        F0 = prctile(C_bg,options.df_prctile,2);
    end
    C_df = bsxfun(@times,Df,1./F0(:));
else
    if options.df_prctile == 50
        if verLessThan('matlab','2015b')
            warning('Median filtering at the boundaries might be inaccurate due to zero padding.')
            F0 = medfilt1(C_bg,options.df_window,[],2);
        else
            F0 = medfilt1(C_bg,options.df_window,[],2,'truncate');
        end
    else
        F0 = zeros(size(C_bg));
        for i = 1:size(Df,1);
            df_temp = running_percentile(C_bg(i,:), options.df_window, options.df_prctile);
            F0(i,:) = df_temp(:)';
        end
    end
    C_df = Df./F0;
end           