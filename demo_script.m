clear;
%% load file

addpath(genpath('../constrained-foopsi'));
addpath(genpath('utilities'));

nam = 'demoMovie.tif';             
                                % insert path to tiff stack here
sframe=1;						% user input: first frame to read (optional, default 1)
num2read=2000;					% user input: how many frames to read   (optional, default until the end)

Y = bigread2(nam,sframe,num2read);
if ~isa(Y,'double');    Y = double(Y);  end         % convert to double

[d1,d2,T] = size(Y);                                % dimensions of dataset
d = d1*d2;                                          % total number of pixels

%% Interpolate any missing data (just for pre-processing)

if any(isempty(Y))
    Y_interp = interp_missing_data(Y);      % interpolate missing data
    mis_data = find(Y_interp);
    Y(mis_data) = Y_interp(mis_data);       % introduce interpolated values for initialization
else
    Y_interp = sparse(d,T);
end

%% fast initialization of spatial components using the greedyROI
k_sub = 1;                                          % temporal interleaving factor (for possible memory issues)

nr = 30;                                           % number of components to be found
params.gSiz = 9;                                   % maximum size of spatial footprint (box of size param.gSiz x param.gSiz)
params.gSig = 4;                                   % std of gaussian (size of neuron)        
[Ain, Cin, center, ~] = greedyROI2d(Y(:,:,1:k_sub:end), nr, params);  
Ain = sparse(reshape(Ain,d,nr));                    % initial estimate of spatial footprints (size d x nr)
Cin = Cin';                                         % 

% display centers of found components
Cn =  correlation_image(Y); %max(Y,[],3); %std(Y,[],3); % image statistic (only for display purposes)
figure;imagesc(Cn);
    axis equal; axis tight; hold all;
    scatter(center(:,2),center(:,1),'mo');
    title('Center of ROIs found from initialization algorithm');
    drawnow;
%% compute estimates of noise for every pixel and a global time constant

Yr = reshape(Y,d,T);
clear Y;
p = 2;                                                          % order of autoregressive system (p=1 just decay, p = 2, both rise and decay)
active_pixels = find(sum(Ain,2));                               % pixels where the greedy method found activity
unsaturated_pixels = find_unsaturatedPixels(Yr);                % pixels that do not exhibit saturation
options.pixels = intersect(active_pixels,unsaturated_pixels);   % base estimates only on                 

P = arpfit(Yr,p,options);
[bin,fin] = nnmf(max(Yr-Ain*Cin,0),1);
P.interp = Y_interp;
P.unsaturatedPix = unsaturated_pixels;
% remove interpolated values
miss_data_int = find(Y_interp);
Yr(miss_data_int) = P.interp(miss_data_int);

%% update spatial components
P.search_method = 'ellipse';
P.d1 = d1;                  % dimensions of image
P.d2 = d2; 
P.dist = 3;                 % ellipse expansion factor for local search of spatial components
[A,b] = update_spatial_components(Yr,Cin,fin,Ain,P);

%% update temporal components
P.method = 'constrained_foopsi';            % choice of method for deconvolution
P.temporal_iter = 2;                        % number of iterations for block coordinate descent
P.fudge_factor = 0.98;                      % fudge factor to reduce time constant estimation bias
[C,f,Y_res,P] = update_temporal_components(Yr,A,b,Cin,fin,P);

%% merge found components

P.merge_thr = 0.8;                          % merging threshold
[Am,Cm,nr_m,merged_ROIs,P] = merge_ROIs(Y_res,A,b,C,f,P);

display_merging = 1; % flag for displaying merging example
if display_merging
    i = 1; randi(length(merged_ROIs));
    ln = length(merged_ROIs{i});
    figure;
        set(gcf,'Position',[300,300,(ln+2)*300,300]);
        for j = 1:ln
            subplot(1,ln+2,j); imagesc(reshape(A(:,merged_ROIs{i}(j)),d1,d2)); 
                title(sprintf('Component %i',j),'fontsize',16,'fontweight','bold'); axis equal; axis tight;
        end
        subplot(1,ln+2,ln+1); imagesc(reshape(Am(:,nr_m-length(merged_ROIs)+i),d1,d2));
                title('Merged Component','fontsize',16,'fontweight','bold');axis equal; axis tight; 
        subplot(1,ln+2,ln+2);
            plot(1:T,(diag(max(C(merged_ROIs{i},:),[],2))\C(merged_ROIs{i},:))'); 
            hold all; plot(1:T,Cm(nr_m-length(merged_ROIs)+i,:)/max(Cm(nr_m-length(merged_ROIs)+i,:)),'--k')
            title('Temporal Components','fontsize',16,'fontweight','bold')
        drawnow;
end

%% repeat
[A2,b2] = update_spatial_components(Yr,Cm,f,Am,P);
[C2,f2,Y_res,P] = update_temporal_components(Yr,A2,b2,Cm,f,P);
C_df = extract_DF_F(Yr,[A2,b2],[C2;f2],nr_m+1); % extract DF/F values (optional)

%% do some plotting

[A_or,C_or,P] = order_ROIs(A2,C2,P);      % order components
contour_threshold = 0.95;                 % amount of energy used for each component to construct contour plot
figure;
[Coor,json_file] = plot_contours(A_or,reshape(P.sn,d1,d2),contour_threshold,1); % contour plot of spatial footprints
pause; 
%savejson('jmesh',json_file,'filename');        % optional save json file with component coordinates (requires matlab json library)
view_patches(Yr,A_or,C_or,b2,f2,d1,d2);                                         % display all components

%% make movie

param.skip_frame = 2;
param.ind = [1,2,3,5];
param.sx = 16;
param.make_avi = 0;
make_patch_video(A_or,C_or,b2,f2,Yr,d1,d2,param)
