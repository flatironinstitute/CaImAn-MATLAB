function upSample(obj, obj_ds, T)
%% upsample the CNMF results 
% input:
%   obj_ds:   class @Sources2D, results on downsampled data 
%   T:        original number of frames 

%% Author: Pengcheng Zhou, Carnegie Mellon University, 2015 

% dimension of raw data
d1 = obj.options.d1;
d2 = obj.options.d2;

% dimension of downsampled data
d1s = obj_ds.options.d1;
d2s = obj_ds.options.d2;
Ts = size(obj_ds.C, 2);

% upsample spatial components
As = reshape(full(obj_ds.A), d1s, d2s, []);
obj.A = reshape(imresize(As, [d1, d2]), d1*d2, []);

bs = reshape(obj_ds.b, d1s, d2s, []);
obj.b = reshape(imresize(bs, [d1, d2]), d1*d2, []);

% update temporal components
if T~=Ts
    Cs = obj_ds.C;
    obj.C = imresize(Cs, [size(Cs, 1), T]);
    
    f = obj_ds.f;
    obj.f = imresize(f, [size(f, 1), T]);
end
% upsample contour coordinate

% update Df

% update C_df

% update S_df

% update options