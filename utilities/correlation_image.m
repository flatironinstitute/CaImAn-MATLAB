function Cn = correlation_image(Y,sz,d1,d2)

% construct correlation image based on neighboing pixels
% Y: raw data
% sz: size of neighborhood (sz either 4 or 8, default 4)
% d1,d2: spatial dimensions

% Author: Eftychios A. Pnevmatikakis, Simons Foundation, 2015

if nargin == 1 || isempty(sz)
    sz = 4;
end

if ismatrix(Y)
    Y = reshape(Y,d1,d2,size(Y,2));
end

[d1,d2,T] = size(Y);

mY = mean(Y,3);
sY = std(Y,[],3);

Y = bsxfun(@minus,Y,mY);

Ycor = zeros(d1,d2,sz);
MASK = zeros(d1,d2);
Yx = sum(Y(1:end-1,:,:).*Y(2:end,:,:),3)./((T-1).*sY(1:end-1,:).*sY(2:end,:));
Ycor(1:end-1,:,1) = Yx;
MASK(1:end-1,:) = MASK(1:end-1,:) + 1;
Ycor(2:end,:,2) = Yx;
MASK(2:end,:) = MASK(2:end,:) + 1;
Yy = sum(Y(:,1:end-1,:).*Y(:,2:end,:),3)./((T-1).*sY(:,1:end-1).*sY(:,2:end));
Ycor(:,1:end-1,3) = Yy;
MASK(:,1:end-1) = MASK(:,1:end-1) + 1;
Ycor(:,2:end,4) = Yy;
MASK(:,2:end) = MASK(:,2:end) + 1;

if sz == 8
    Yxpyp = sum(Y(1:end-1,1:end-1,:).*Y(2:end,2:end,:),3)./((T-1).*sY(1:end-1,1:end-1).*sY(2:end,2:end));
    Ycor(1:end-1,1:end-1,5) = Yxpyp;
    MASK(1:end-1,1:end-1) = MASK(1:end-1,1:end-1) + 1;
    Ycor(2:end,2:end,6) = Yxpyp;
    MASK(2:end,2:end) = MASK(2:end,2:end) + 1;
    Yxpym = sum(Y(1:end-1,2:end,:).*Y(2:end,1:end-1,:),3)./((T-1).*sY(1:end-1,2:end).*sY(2:end,1:end-1));
    Ycor(1:end-1,2:end,7) = Yxpym;
    MASK(1:end-1,2:end) = MASK(1:end-1,2:end) + 1;
    Ycor(2:end,1:end-1,8) = Yxpym;    
    MASK(2:end,1:end-1) = MASK(2:end,1:end-1) + 1;
end

Cn = sum(Ycor,3)./MASK;