function [A,b,C] = update_spatial_components(Y,C,f,A_,P,options)

% update spatial footprints and background through Basis Pursuit Denoising
% for each pixel i solve the problem 
%   [A(i,:),b(i)] = argmin sum(A(i,:))
%       subject to || Y(i,:) - A(i,:)*C + b(i)*f || <= sn(i)*sqrt(T);
% for each pixel the search is limited to a few spatial components

% INPUTS:
% Y:    raw data
% C:    temporal components
% f:    temporal background
% A_:   current estimate of spatial footprints (used for determining search locations only) 
% P:    dataset parameters (used for noise values and interpolated entries)

% options    parameter struct (for noise values and other parameters)

% OUTPUTS:
% A:    new estimate of spatial footprints
% b:    new estimate of spatial background
% C:    temporal components (updated only when spatial components are completely removed)

% Written by:
% Eftychios A. Pnevmatikakis, Simons Foundation, 2015

warning('off', 'MATLAB:maxNumCompThreads:Deprecated');
memmaped = isobject(Y);
if memmaped
    sizY = size(Y,'Y');
    d = prod(sizY(1:end-1));
    T = sizY(end);
else
    [d,T] = size(Y);
end
if nargin < 6 || isempty(options); options = []; end
if ~isfield(options,'d1') || isempty(options.d1); d1 = input('What is the total number of rows? \n'); options.d1 = d1; else d1 = options.d1; end          % # of rows
if ~isfield(options,'d2') || isempty(options.d2); d2 = input('What is the total number of columns? \n'); options.d2 = d2; else d2 = options.d2; end          % # of columns
if ~isfield(options,'show_sum'); show_sum = 0; else show_sum = options.show_sum; end            % do some plotting while calculating footprints
if ~isfield(options,'interp'); Y_interp = sparse(d,T); else Y_interp = options.interp; end      % identify missing data
if ~isfield(options,'use_parallel'); use_parallel = ~isempty(which('parpool')); else use_parallel = options.use_parallel; end % use parallel toolbox if present
if ~isfield(options,'search_method'); method = []; else method = options.search_method; end     % search method for determining footprint of spatial components

if nargin < 2 || (isempty(A_) && isempty(C))  % at least either spatial or temporal components should be provided
    error('Not enough input arguments')
else
    if ~isempty(C); K = size(C,1); elseif islogical(A_); K = size(A_,2); else K = size(A_2,2) - options.nb; end
end

if nargin < 5 || isempty(P); P = preprocess_data(Y,1); end  % etsimate noise values if not present
if nargin < 4 || isempty(A_); 
    IND = ones(d,size(C,1)); 
else
    if islogical(A_)     % check if search locations have been provided, otherwise estimate them
        IND = A_;
        if isempty(C)    
            INDav = double(IND)/diag(sum(double(IND)));          
            px = (sum(IND,2)>0);
            f = mean(Y(~px,:));
            b = max(Y*f',0)/norm(f)^2;
            C = max(INDav'*Y - (INDav'*b)*f,0);
        end
    else
        IND = determine_search_location(A_(:,1:K),method,options);
    end
end

options.sn = P.sn;
if ~memmaped
    Y(P.mis_entries) = NaN; % remove interpolated values
end

Cf = [C;f];

if use_parallel         % solve BPDN problem for each pixel
    Nthr = max(2*maxNumCompThreads,round(d*T/2^24));
    siz_row = [floor(d/Nthr)*ones(Nthr-1,1);d-floor(d/Nthr)*(Nthr-1)];
    indeces = [0;cumsum(siz_row)];
    if ~memmaped
        Ycell = mat2cell(Y,siz_row,T);
    else
        Ycell = cell(Nthr,1);
    end
    INDc =  mat2cell(IND,siz_row,K);
    Acell = cell(Nthr,1);
    Psnc = mat2cell(options.sn(:),siz_row,1); 
    Yf = cell(Nthr,1);
    parfor nthr = 1:Nthr
        Acell{nthr} = zeros(siz_row(nthr),size(Cf,1));
        if memmaped
            Ytemp = double(Y.Yr(indeces(nthr)+1:indeces(nthr+1),:));
        else
            Ytemp = Ycell{nthr};
        end
        Yf{nthr} = Ytemp*f';
        for px = 1:siz_row(nthr)
            fn = ~isnan(Ytemp(px,:));       % identify missing data
            ind = find(INDc{nthr}(px,:));
            if ~isempty(ind);
                ind2 = [ind,K+(1:size(f,1))];
                %[~, ~, a, ~] = lars_regression_noise(Ycell{nthr}(px,fn)', Cf(ind2,fn)', 1, Psnc{nthr}(px)^2*T);
                [~, ~, a, ~] = lars_regression_noise(Ytemp(px,fn)', Cf(ind2,fn)', 1, Psnc{nthr}(px)^2*T);
                Acell{nthr}(px,ind2) = a';
            end
        end
    end
    A = cell2mat(Acell);
    Yf = cell2mat(Yf);
else
    A = [zeros(d,K),zeros(d,size(f,1))];
    sA = zeros(d1,d2);
    Yf = Y*f';
    for px = 1:d   % estimate spatial components
        fn = ~isnan(Y(px,:));       % identify missing data
        ind = find(IND(px,:));
        if ~isempty(ind);
            ind2 = [ind,K+(1:size(f,1))];
            [~, ~, a, ~] = lars_regression_noise(Y(px,fn)', Cf(ind2,fn)', 1, options.sn(px)^2*T);
            A(px,ind2) = a';
            sA(px) = sum(a);
        end
        if show_sum
            if mod(px,d1) == 0;
               figure(20); imagesc(sA); axis square;  
               title(sprintf('Sum of spatial components (%i out of %i columns done)',round(px/d1),d2)); drawnow;
            end
        end
    end
end

A(isnan(A))=0;
A = sparse(A);
A = threshold_components(A,options);  % post-processing of components

fprintf('Updated spatial components \n');

ff = find(sum(A(:,1:K))==0);           % remove empty components
if ~isempty(ff)
    K = K - length(ff);
    A(:,ff) = [];
    C(ff,:) = [];
end

b = max((Yf - A(:,1:K)*(C(1:K,:)*f'))/(f*f'),0);
A = A(:,1:K);