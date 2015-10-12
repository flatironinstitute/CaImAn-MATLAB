function [A,b] = update_spatial_components(Y,C,f,A_,P)

% update spatial footprints and background through Basis Pursuit Denoising
% for each pixel i solve the problem 
%   [A(i,:),b(i)] = argmin sum(A(i,:))
%       subject to || Y(i,:) - A(i,:)*C + b(i)*f || <= P.sn(i)*sqrt(T);
% for each pixel the search is limited to a few spatial components

% INPUTS:
% Y:    raw data
% C:    temporal components
% f:    temporal background
% A_:   current estimate of spatial footprints (used for determining search locations only)
% P:    parameter vector (for noise values and other parameters)

% OUTPUTS:
% A:    new estimate of spatial footprints
% b:    new estimate of spatial background

% Written by:
% Eftychios A. Pnevmatikakis, Simons Foundation, 2015

warning('off', 'MATLAB:maxNumCompThreads:Deprecated');

[d,T] = size(Y);
if ~isfield(P,'d1'); d1 = sqrt(d); else d1 = P.d1; end          % # of rows
if ~isfield(P,'d2'); d2 = sqrt(d); else d2 = P.d2; end          % # of columns
if ~isfield(P,'min_size'); min_size = 8; else min_size = P.min_size; end    % minimum size of ellipse axis
if ~isfield(P,'max_size'); max_size = 3; else max_size = P.max_size; end    % maximum size of ellipse axis
if ~isfield(P,'dist'); dist = 3; else dist = P.dist; end                    % expansion factor of ellipse
if ~isfield(P,'show_sum'); show_sum = 0; else show_sum = P.show_sum; end            % do some plotting while calculating footprints
if ~isfield(P,'interp'); Y_interp = sparse(d,T); else Y_interp = P.interp; end      % identify missing data
if ~isfield(P,'use_parallel'); use_parallel = ~isempty(which('parpool')); else use_parallel = P.use_parallel; end % use parallel toolbox if present
if ~isfield(P,'search_method'); method = []; else method = P.search_method; end     % search method for determining footprint of spatial components

nr = size(C,1);       % number of neurons

IND = determine_search_location(A_(:,1:nr),method,P);

Cf = [C;f];

if use_parallel         % solve BPDN problem for each pixel
    Nthr = 2*maxNumCompThreads;
    siz_row = [ceil(d/Nthr)*ones(Nthr-1,1);d-ceil(d/Nthr)*(Nthr-1)];    
    Ycell = mat2cell(Y,siz_row,T);
    INDc =  mat2cell(IND,siz_row,nr);
    Acell = cell(Nthr,1);
    Psnc = mat2cell(P.sn,siz_row,1);    
    parfor nthr = 1:Nthr
        Acell{nthr} = zeros(siz_row(nthr),size(Cf,1));
        for px = 1:siz_row(nthr)
            fn = ~isnan(Ycell{nthr}(px,:));       % identify missing data
            ind = find(INDc{nthr}(px,:));
            if ~isempty(ind);
                ind2 = [ind,nr+(1:size(f,1))];
                [~, ~, a, ~] = lars_regression_noise(Ycell{nthr}(px,fn)', Cf(ind2,fn)', 1, Psnc{nthr}(px)^2*T);
                Acell{nthr}(px,ind2) = a';
            end
        end
    end
    A = cell2mat(Acell);
else
    A = [zeros(d,nr),zeros(d,size(f,1))];
    sA = zeros(d1,d2);
    for px = 1:d   % estimate spatial components
        fn = ~isnan(Y(px,:));       % identify missing data
        ind = find(IND(px,:));
        if ~isempty(ind);
            ind2 = [ind,nr+(1:size(f,1))];
            [~, ~, a, ~] = lars_regression_noise(Y(px,fn)', Cf(ind2,fn)', 1, P.sn(px)^2*T);
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
A = threshold_components(A,P);  % post-processing of components

fprintf('Updated spatial components \n');

% ff = find(sum(A)==0);           % remove empty components
% if ~isempty(ff)
%     nr = nr - length(ff);
%     A(:,ff) = [];
%     C(ff,:) = [];
% end

if nnz(Y_interp);
    ff = find(Y_interp);
    Y(ff) = Y_interp(ff);
end

Y_res = Y - A(:,1:nr)*C(1:nr,:);
A_bas = max(Y_res*f'/norm(f)^2,0); % update baseline based on residual
b = A_bas;
A = A(:,1:nr);