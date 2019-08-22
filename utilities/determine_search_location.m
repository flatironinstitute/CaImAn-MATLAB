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
    method = 'dilate';         % default method
end

[d,nr] = size(A);

if ~isfield(params,'d1'); d1 = sqrt(d); params.d1 = d1; else d1 = params.d1; end          % # of rows
if ~isfield(params,'d2'); d2 = sqrt(d); params.d2 = d2; else d2 = params.d2; end          % # of columns
if ~isfield(params,'d3'); d3 = 1; params.d3 = d3; else d3 = params.d3; end
if ~isfield(params,'min_size') || isempty(params.min_size); min_size = 3; else min_size = params.min_size; end    % minimum size of ellipse axis 
if ~isfield(params,'max_size') || isempty(params.max_size); max_size = 8; else max_size = params.max_size; end    % maximum size of ellipse axis
if ~isfield(params,'dist')  || isempty(params.dist); dist = 3; else dist = params.dist; end                              % expansion factor of ellipse
if ~isfield(params,'se')  || isempty(params.se);  % morphological element (for 'dilate')
    if d3 == 1; 
        expandCore = strel('disk',4,0); 
    else
        [xx,yy,zz] = meshgrid(-2:2,-2:2,-1:1);
        expandCore = strel(xx.^2 + yy.^2 + 2*zz.^2 <= 2^2);
    end
else    
    expandCore = params.se;     
    if d3 > 1 && ismatrix(expandCore.getnhood)
        [xx,yy,zz] = meshgrid(-2:2,-2:2,-1:1);
        expandCore = strel(xx.^2 + yy.^2 + 2*zz.^2 <= 2^2);
    end
end     

if strcmpi(method,'ellipse'); method = 'ellipse';
elseif strcmpi(method,'dilate'); method = 'dilate';
else fprintf('Method not recongnized. Search location equals the entire field of view. \n');
end

IND = logical(sparse(d,nr));
if strcmpi(method,'ellipse') && (dist == Inf)
   IND = true(d,nr);
else
    cm = com(A,d1,d2,d3);
    if d3 == 1
        cm = [cm,ones(nr,1)];
    end
    %if strcmpi(method,'ellipse'); Vr = cell(nr,1); end
    parfor i = 1:nr            % calculation of variance for each component and construction of ellipses
        sx = 32;
        ind_x = round([max(1,cm(i,1)-sx),min(d1,cm(i,1)+sx)]);
        ind_y = round([max(1,cm(i,2)-sx),min(d2,cm(i,2)+sx)]);
        ind_z = round([max(1,cm(i,3)-sx),min(d3,cm(i,3)+sx)]);
        siz_p = [ind_x(2)-ind_x(1)+1,ind_y(2)-ind_y(1)+1,ind_z(2)-ind_z(1)+1]; 
        dp = prod(siz_p);
       idx = patch_to_linear([ind_x,ind_y,ind_z], [d1,d2,d3,1]);
       switch method
           case 'ellipse'
               cx = kron(ones(siz_p(2)*siz_p(3),1),(1:siz_p(1))'); 
               cy = kron(ones(siz_p(3),1),kron((1:siz_p(2))',ones(siz_p(1),1)));
               cz = kron((1:siz_p(3))',ones(siz_p(1)*siz_p(2),1));
               if d3 == 1
                   Vr = ([cx - (cm(i,1) - ind_x(1) + 1), cy - (cm(i,2) - ind_y(1) + 1)]'*spdiags(A(idx,i),0,dp,dp)*[cx - (cm(i,1) - ind_x(1) + 1), cy - (cm(i,2) - ind_y(1) + 1)])/sum(A(idx,i));
                   [V,D] = eig(Vr);
                   cor = [cx - (cm(i,1) - ind_x(1) + 1), cy - (cm(i,2) - ind_y(1) + 1)];
               else
                   Vr = ([cx - (cm(i,1) - ind_x(1) + 1), cy - (cm(i,2) - ind_y(1) + 1), cz - (cm(i,3) - ind_z(1) + 1)]'*spdiags(A(idx,i),0,dp,dp)*[cx - (cm(i,1) - ind_x(1) + 1), cy - (cm(i,2) - ind_y(1) + 1), cz - (cm(i,3) - ind_z(1) + 1)])/sum(A(idx,i));
                   [V,D] = eig(Vr);
                   cor = [cx - (cm(i,1) - ind_x(1) + 1), cy - (cm(i,2) - ind_y(1) + 1), cz - (cm(i,3) - ind_z(1) + 1)];
               end                              
               d11 = min(max_size^2,max(min_size^2,D(1,1)));
               d22 = min(max_size^2,max(min_size^2,D(2,2)));               
               if d3 == 1
                   ind_t = sqrt((cor*V(:,1)).^2/d11 + (cor*V(:,2)).^2/d22)<=dist;
               else
                   d33 = min((max_size/2)^2,max((min_size/2)^2,D(3,3)));
                   ind_t = sqrt((cor*V(:,1)).^2/d11 + (cor*V(:,2)).^2/d22 + (cor*V(:,3)).^2/d33)<=dist;       % search indeces for each component
               end
               ind_temp = sparse(idx,1,ind_t,d,1); 
               IND(:,i) = logical(ind_temp); 
           case 'dilate'
               A_temp = imdilate(reshape(full(A(:,i)), [d1,d2,d3]), expandCore);
               IND(:,i) = A_temp(:)>0;
       end
   end
end


end
