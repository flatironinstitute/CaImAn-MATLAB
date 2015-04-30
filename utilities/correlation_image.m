function Cn = correlation_image(Y,sz,d1,d2)

% construct correlation image based on neighboing pixels
% Y: raw data
% sz: size of neighborhood (sz either 4 or 8, default 4)
% d1,d2: spatial dimensions

if nargin == 1 || isempty(sz)
    sz = 4;
end

if ndims(Y) == 3
    [d1,d2,~] = size(Y);
    Y = reshape(Y,d1*d2,size(Y,3));
end

d = d1*d2;
Cn = zeros(d1,d2);

if sz == 8
    NB = [ 1 0; -1 0; 1 1; 0 1; -1 1; 1 -1; 0 -1; -1 -1];
else
    NB = [ 1 0; -1 0; 0 1; 0 -1];
end

if sz == 4
    if mod(d1,2)
        A = false(d1,d2);
        A(1:2:end) = 1;
    else
        A = false(d1+1,d2);
        A(1:2:end) = 1;
        A = A(1:d1,:);
    end
else
    temp = zeros(d1,1);
    temp(2:2:end) = 1;
    A = repmat([ones(d1,1),temp],1,floor(d2/2));
    if mod(d2,2)
        A = [A,ones(d1,1)];
    end
end

ff = find(A);
%CM = cell(length(ff),1);
CM = cell(d,1);
for k = 1:length(ff)
    i = ff(k);
    loc = [i - d1*(ceil(i/d1)-1), ceil(i/d1)];
    nb_loc = NB + repmat(loc,[sz 1]);
        nb_loc(nb_loc(:,1)>d1,:) = [];
        nb_loc(nb_loc(:,1)<1,:) = [];
        nb_loc(nb_loc(:,2)>d2,:) = [];
        nb_loc(nb_loc(:,2)<1,:) = [];
        nb_loc = nb_loc(:,1) + d1*(nb_loc(:,2)-1);
    temp_corr = corr(Y(i,:)',Y(nb_loc,:)');
    %CM{k} = [nb_loc(:),temp_corr(:)];
    CM{i} = [nb_loc(:),temp_corr(:)];
    Cn(loc(1),loc(2)) = mean(temp_corr);
end

cc = setdiff(1:d,ff);
for k = 1:length(cc)
    i = cc(k);
    loc = [i - d1*(ceil(i/d1)-1), ceil(i/d1)];
    nb_loc = nb_location(loc);
    temp = zeros(size(nb_loc));
    for j = 1:length(nb_loc)
        %temp(j) = CM{ff==nb_loc(j)}(CM{ff==nb_loc(j)}(:,1)==i,2);
        temp(j) = CM{nb_loc(j)}(CM{nb_loc(j)}(:,1)==i,2);
    end
    Cn(loc(1),loc(2)) = mean(temp);
end

    function nb_loc = nb_location(loc)
        nb_loc = NB + repmat(loc,[sz 1]);
        nb_loc(nb_loc(:,1)>d1,:) = [];
        nb_loc(nb_loc(:,1)<1,:) = [];
        nb_loc(nb_loc(:,2)>d2,:) = [];
        nb_loc(nb_loc(:,2)<1,:) = [];
        nb_loc = nb_loc(:,1) + d1*(nb_loc(:,2)-1);
    end

end