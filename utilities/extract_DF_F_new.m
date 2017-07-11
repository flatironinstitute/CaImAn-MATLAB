function [C_df,Df] = extract_DF_F_new(A,C,b,f,P,options)

% extract DF/F signals after performing NMF
% inputs:  Y raw data (d X T matrix, d # number of pixels, T # of timesteps)
%          A matrix of spatial components (d x K matrix, K # of components)
%          C matrix of temporal components (K x T matrix)
%          P neuron structure, used to read the baseline activity for each
%                    component of C
%          options structure used for specifying method for determining DF
%           default method is the median of the trace. By changing
%           options.df_prctile an arbitray percentile can be used (between 0 and 100).
%           a moving window can also be established by specifying options.df_window

% outputs:  C_df temporal components in the DF/F domain
%           Df   background for each component to normalize the filtered raw data    

% Written by: 
% Eftychios A. Pnevmatikakis, Simons Foundation, 2016

%memmaped = isobject(Y);
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
if ~isfield(options,'full_A') || isempty(options.full_A); full_A = defoptions.full_A; else full_A = options.full_A; end

[K,T] = size(C);

Bas = zeros(K,T);

if ~(nargin < 5 || isempty(P))
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

nA = sqrt(sum(A.^2));

AA = A'*A;
AA(1:K+1:end) = 0;

Cf = bsxfun(@times,C - Bas,nA(:).^2);
C2 = repmat(AA*bas_val,1,T) + (A'*b)*f;

if isempty(options.df_window) || (options.df_window > size(C,2))
    if options.df_prctile == 50
        Df = median(C2,2);
    else
        Df = prctile(C2,options.df_prctile,2);
    end
    C_df = bsxfun(@times,Cf,1./Df(:));
else
    if options.df_prctile == 50
        if verLessThan('matlab','2015b')
            warning('Median filtering at the boundaries might be inaccurate due to zero padding.')
            Df = medfilt1(C2,options.df_window,[],2);
        else
            Df = medfilt1(C2,options.df_window,[],2,'truncate');
        end
    else
        Df = zeros(size(C2));
        for i = 1:size(Df,1);
            df_temp = running_percentile(C2(i,:), options.df_window, options.df_prctile);
            Df(i,:) = df_temp(:)';
        end
    end
    C_df = Cf./Df;
end           