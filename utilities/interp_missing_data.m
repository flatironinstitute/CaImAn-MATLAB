function Y_interp = interp_missing_data(Y)

% interpolate missing data using linear interpolation for each pixel
% produce a sparse matrix with the values 

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