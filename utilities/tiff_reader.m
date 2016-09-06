function Y = tiff_reader(name,T)

% read tiff stack. Optional truncate to first T timesteps 

info = imfinfo(name);

if nargin == 1
   T = numel(info);
end
    
d1 = info(1).Height;
d2 = info(1).Width;

Y = zeros(d1,d2,T);
for t = 1:T
    Y(:,:,t) = imread(name, t, 'Info',info);
end