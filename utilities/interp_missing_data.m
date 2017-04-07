function Y_interp = interp_missing_data(Y)
% INTERP_MISSING_DATA - interpolate missing data using linear interpolation for each pixel
%
%  [Y_INTERP] = INTERP_MISSING_DATA(Y)
%
%  Given a matrix Y with possible NaN values, this function produces
%  a sparse matrix Y_INTERP that has the linearly interpolated values
%  of each pixel that exhibits a NaN value. The values are interpolated
%  over the last dimension of Y. 
%
%  Example:
%     yy = [0 0 ; 0 0];
%     yy(:,:,2) = [0 NaN ; 0 2];
%     yy(:,:,3) = [0.1 0.2 ; 0 2.3];
%     Y_interp = interp_missing_data(yy);
%     % Y_interp is a sparse matrix :    (3,2)  0.1000
%
%  See also: INTERP1

sizY = size(Y);
dimY = length(sizY);
d = prod(sizY(1:dimY-1));
T = sizY(end);
mis_data = cell(d,1);

for i = 1:d
    [ii,jj,kk] = ind2sub(sizY(1:dimY-1),i);
    if dimY == 2
        ytemp = Y(i,:);
    elseif dimY == 3
        ytemp = squeeze(Y(ii,jj,:));
    elseif dimY == 4
        ytemp = squeeze(Y(ii,jj,kk,:));
    end
    f = isnan(ytemp(:));
    y_val = interp1(find(~f),ytemp(~f),find(f),'linear','extrap');
    mis_data{i} = [i*ones(length(y_val),1),find(f(:)),y_val(:)];
end

mis_data = cell2mat(mis_data);

Y_interp = sparse(mis_data(:,1),mis_data(:,2),mis_data(:,3),d,T);
