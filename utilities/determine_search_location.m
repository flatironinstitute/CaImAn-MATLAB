function IND = determine_search_location(A,method,params)

% determine the search location for updating each spatial component
% INPUTS:
% A:        current estimate of spatial component ( d x nr sparse matrix)
% method:   method to be used of determining search locations
%     available methods:
%       'ellipse': draw an ellipse centered at the center of mass with
%           axes the first two principal components and expand it to
%           determine the search location  (default)
%       'dilate' : grow the spatial footprint by dilating the current
%           footprint
%      if a different argument is passed then search location covers the whole field of view
% params:   hyper-parameter struct for the various methods (see default settings below for description)
%
% OUTPUT:
% IND:      binary d x nr matrix. IND(i,j) = 1 if pixel i is included in the search location of component j
%
% Written by Eftychios A. Pnevmatikakis, Simons Foundation 
%   with input from Weijian Yang, Columbia University

if nargin < 3
    params = [];
end

if nargin < 2 || isempty(method)
    method = 'ellipse';         % default method
end

[d,nr] = size(A);
if ~isfield(params,'d1'); d1 = sqrt(d); params.d1 = d1; else d1 = params.d1; end          % # of rows
if ~isfield(params,'d2'); d2 = sqrt(d); params.d2 = d2; else d2 = params.d2; end          % # of columns
if ~isfield(params,'min_size') || isempty(params.min_size); min_size = 3; else min_size = params.min_size; end    % minimum size of ellipse axis 
if ~isfield(params,'max_size') || isempty(params.max_size); max_size = 8; else max_size = params.max_size; end    % maximum size of ellipse axis
if ~isfield(params,'dist')  || isempty(params.dist); dist = 3; else dist = params.dist; end                              % expansion factor of ellipse
if ~isfield(params,'se')  || isempty(params.se); expandCore = strel('disk',4,0); else expandCore = params.se; end      % morphological element (for 'dilate')

if strcmpi(method,'ellipse'); method = 'ellipse';
elseif strcmpi(method,'dilate'); method = 'dilate';
else fprintf('Method not recongnized. Search location equals the entire field of view. \n');
end

IND = false(d,nr);
switch method 
    case 'ellipse'
        Coor.x = kron(ones(d2,1),(1:d1)'); 
        Coor.y = kron((1:d2)',ones(d1,1));
        if ~(dist==Inf)             % determine search area for each neuron
           cm = zeros(nr,2);        % vector for center of mass
           Vr = cell(nr,1);
           IND = zeros(d,nr);       % indicator for distance								   
           cm(:,1) = Coor.x'*A(:,1:nr)./sum(A(:,1:nr)); 
           cm(:,2) = Coor.y'*A(:,1:nr)./sum(A(:,1:nr));          % center of mass for each components
           for i = 1:nr            % calculation of variance for each component and construction of ellipses
               Vr{i} = ([Coor.x - cm(i,1), Coor.y - cm(i,2)]'*spdiags(A(:,i),0,d,d)*[Coor.x - cm(i,1), Coor.y - cm(i,2)])/sum(A(:,i));
               [V,D] = eig(Vr{i});
               cor = [Coor.x - cm(i,1),Coor.y - cm(i,2)];
               d11 = min(max_size^2,max(min_size^2,D(1,1)));
               d22 = min(max_size^2,max(min_size^2,D(2,2)));
               IND(:,i) = sqrt((cor*V(:,1)).^2/d11 + (cor*V(:,2)).^2/d22)<=dist;       % search indeces for each component
           end
        else
            IND = true(d,nr);
        end
    case 'dilate'
        A = threshold_components(A,params);
        for i = 1:nr
            A_temp = imdilate(reshape(full(A(:,i)),d1,d2),expandCore);
            IND(:,i) = A_temp(:)>0;
        end
    otherwise
        IND = true(d,nr);
end
