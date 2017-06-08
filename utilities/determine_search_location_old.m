function IND = determine_search_location_old(A,method,params)

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
if ~isfield(params,'d3'); d3 = 1; params.d3 = d3; else d3 = params.d3; end
if ~isfield(params,'min_size') || isempty(params.min_size); min_size = 3; else min_size = params.min_size; end    % minimum size of ellipse axis 
if ~isfield(params,'max_size') || isempty(params.max_size); max_size = 8; else max_size = params.max_size; end    % maximum size of ellipse axis
if ~isfield(params,'dist')  || isempty(params.dist); dist = 3; else dist = params.dist; end                              % expansion factor of ellipse
if ~isfield(params,'se')  || isempty(params.se);  % morphological element (for 'dilate')
    if d3 == 1; expandCore = strel('disk',4,0); else expandCore = strel(ones(4,4,2)); end
else
    expandCore = params.se; 
end     

if strcmpi(method,'ellipse'); method = 'ellipse';
elseif strcmpi(method,'dilate'); method = 'dilate';
else fprintf('Method not recongnized. Search location equals the entire field of view. \n');
end

IND = logical(sparse(d,nr));
switch method 
    case 'ellipse'
        Coor.x = kron(ones(d2*d3,1),(1:d1)'); 
        Coor.y = kron(ones(d3,1),kron((1:d2)',ones(d1,1)));
        Coor.z = kron((1:d3)',ones(d1*d2,1));
        if ~(dist==Inf)             % determine search area for each neuron
           %cm = zeros(nr,2);        % vector for center of mass
           cm = com(A,d1,d2,d3);
           if d3 == 1
               cm = [cm,ones(nr,1)];
           end
           Vr = cell(nr,1);
           %cm(:,1) = Coor.x'*A(:,1:nr)./sum(A(:,1:nr)); 
           %cm(:,2) = Coor.y'*A(:,1:nr)./sum(A(:,1:nr));          % center of mass for each components
           parfor i = 1:nr            % calculation of variance for each component and construction of ellipses
               if d3 == 1
                   Vr{i} = ([Coor.x - cm(i,1), Coor.y - cm(i,2)]'*spdiags(A(:,i),0,d,d)*[Coor.x - cm(i,1), Coor.y - cm(i,2)])/sum(A(:,i));
                   [V,D] = eig(Vr{i});
                   cor = [Coor.x - cm(i,1),Coor.y - cm(i,2)];
               else
                   Vr{i} = ([Coor.x - cm(i,1), Coor.y - cm(i,2), Coor.z - cm(i,3)]'*spdiags(A(:,i),0,d,d)*[Coor.x - cm(i,1), Coor.y - cm(i,2), Coor.z - cm(i,3)])/sum(A(:,i));
                   [V,D] = eig(Vr{i});
                   cor = [Coor.x - cm(i,1),Coor.y - cm(i,2),Coor.z - cm(i,3)];
               end                              
               d11 = min(max_size^2,max(min_size^2,D(1,1)));
               d22 = min(max_size^2,max(min_size^2,D(2,2)));               
               if d3 == 1
                   IND(:,i) = sqrt((cor*V(:,1)).^2/d11 + (cor*V(:,2)).^2/d22)<=dist;
               else
                   d33 = min((max_size/2)^2,max((min_size/2)^2,D(3,3)));
                   IND(:,i) = sqrt((cor*V(:,1)).^2/d11 + (cor*V(:,2)).^2/d22 + (cor*V(:,3)).^2/d33)<=dist;       % search indeces for each component
               end
           end
        else
            IND = true(d,nr);
        end
    case 'dilate'
        A = threshold_components(A,params);
        parfor i = 1:nr
            A_temp = imdilate(reshape(full(A(:,i)),d1,d2,d3),expandCore);
            IND(:,i) = A_temp(:)>0;
        end
    otherwise
        IND = true(d,nr);
end
