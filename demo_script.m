clear;
%% load file

addpath(genpath('utilities'));
nam = ''; % insert path to tiff stack here
Y = tiff_reader(nam);
[d1,d2,T] = size(Y);  
d = d1*d2;

Y_interp = interp_missing_data(Y);      % interpolate missing data (just for pre-processing)
mis_data = find(Y_interp);
Y(mis_data) = Y_inter(mis_data);        % introduce interpolated values for initialization
%% fast initialization of spatial components using the greedyROI

int = 1:4000;            % interval to be processed (due to memory issues)   

nr = 350;                % number of components to be found
params.gSiz = 15;        % maximum size of neuron in pixels
params.gSig = 8;         % std of gaussian (size of neuron)        
[basis, Cin, center, ~] = greedyROI2d(Y(:,:,int), nr, params);  % reduce size for memory reasons
Ain = sparse(reshape(basis,d,nr));  
Cin = Cin';
clear basis;
    % display centers of found components
Cn = mean(Y,3); %correlation_image(Y); %max(Y,[],3); % image statistic
figure;imagesc(Cn);
    axis equal; axis tight; hold all;
    scatter(center(:,2),center(:,1),'mo');
    title('Center of ROIs found from initialization algorithm');
    
%% compute estimates of noise for every pixel and a global time constant

ff = find(sum(Ain,2)); % pixels were greedy method found activity 
p = 1;                 % order of AR system
options.pixels = ff;
Yr = reshape(Y(:,:,int),d,length(int));
P = arpfit(Yr,p,options);
[bin,fin] = nnmf(max(Yr-Ain*Cin,0),1);
P.interp = Y_interp(:,int);
% remove interpolated values
miss_data_int = find(Y_interp(:,int));
Yr(mis_data_int) = P.interp(mis_data_int);

%% update spatial components

[A,b] = update_spatial_components(Yr,Cin,fin,Ain,P);

%% update temporal components
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
[C2,f2,Y_res] = update_temporal_components(Yr,A_or,b,C_or,f,P);
%% do some plotting

[A_or,C_or] = order_ROIs(A,C2);
[Coor,json_file] = plot_contours(A_or,reshape(P.sn,d1,d2),[],1);
view_patches(Yr,A_or,C_or,b,f2,d1,d2)
%savejson('jmesh',json_file,'json-005');

%%
param.skip_frame = 3;
param.ind = [5,8,10,11];
param.make_avi = 1;
make_patch_video(A_or,C2,b,f2,Yr,d1,d2,param)
