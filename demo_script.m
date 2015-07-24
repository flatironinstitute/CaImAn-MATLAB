clear;
%% load file

addpath(genpath('utilities'));
nam = '';             % insert path to tiff stack here
sframe=1;						% user input: first frame to read (optional, default 1)
num2read=1000;					% user input: how many frames to read   (optional, default until the end)

Y = bigread2(nam,sframe,num2read);
if ~isa(Y,'double');    Y = double(Y);  end         % convert to double

[d1,d2,T] = size(Y);                    % dimensions of dataset
d = d1*d2;                              % total number of pixels


%% Interpolate missing data (just for pre-processing)

if any(isempty(Y))
    Y_interp = interp_missing_data(Y);      % interpolate missing data
    mis_data = find(Y_interp);
    Y(mis_data) = Y_interp(mis_data);       % introduce interpolated values for initialization
else
    Y_interp = sparse(d,T);
end

%% fast initialization of spatial components using the greedyROI
k_sub = 1;                                          % temporal interleaving factor (for possible memory issues)

nr = 200;                                           % number of components to be found
params.gSiz = 15;                                   % maximum size of spatial footprint (box of size param.gSiz x param.gSiz)
params.gSig = 8;                                    % std of gaussian (size of neuron)        
[Ain, Cin, center, ~] = greedyROI2d(Y(:,:,1:k_sub:end), nr, params);  
Ain = sparse(reshape(Ain,d,nr));                    % initial estimate of spatial footprints (size d x nr)
Cin = Cin';                                         % 

% display centers of found components
Cn =  mean(Y,3); %correlation_image(Y); %max(Y,[],3); % image statistic (only for display purposes)
figure;imagesc(Cn);
    axis equal; axis tight; hold all;
    scatter(center(:,2),center(:,1),'mo');
    title('Center of ROIs found from initialization algorithm');
    
%% compute estimates of noise for every pixel and a global time constant
                          
p = 1;                                              % order of autoregressive system (p=1 just decay, p = 2, both rise and decay)
options.pixels = find(sum(Ain,2));                  % base estimates only on pixels where the greedy method found activity                
Yr = reshape(Y,d,T);
P = arpfit(Yr,p,options);
[bin,fin] = nnmf(max(Yr-Ain*Cin,0),1);
P.interp = Y_interp;
% remove interpolated values
miss_data_int = find(Y_interp(:,int));
Yr(miss_data_int) = P.interp(miss_data_int);

%% update spatial components
P.d1 = d1; P.d2 = d2; P.dist = 3;
[A,b] = update_spatial_components(Yr,Cin,fin,Ain,P);

%% update temporal components
P.method = 'project';
[C,f,Y_res] = update_temporal_components(Yr,A,b,Cin,fin,P);

%% merge found components
A_in = A;
C_in = C;
repeat = 1;
A_ = A;
C_ = C;
P.method = 'project';
P.merge_thr = 0.8;
while repeat
    [A,C,nr,merged_ROIs] = merge_ROIs(Y_res,A_,b,C_,f,P);
    repeat = ~isempty(merged_ROIs);
    disp(nr)
    A_ = A;
    C_ = C;
end

%% repeat
%[A,b] = update_spatial_components(Yr,C,f,A,P);

P.method = 'constrained_foopsi';
[C2,f2,Y_res] = update_temporal_components(Yr,A_,b,C_,f,P);
%% do some plotting

[A_or,C_or] = order_ROIs(A,C2);
[Coor,json_file] = plot_contours(A_or,reshape(P.sn,d1,d2),[],1);
view_patches(Yr,A_or,C_or,b,f2,d1,d2)
%savejson('jmesh',json_file,'json-005');

%%
param.skip_frame = 3;
param.ind = [1:4];
param.make_avi = 0;
make_patch_video(A_or,C_or,b,f2,Yr,d1,d2,param)
