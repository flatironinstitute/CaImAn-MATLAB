clear;
%% load file

addpath(genpath('../constrained-foopsi'));
addpath(genpath('utilities'));
             
nam = 'demoMovie.tif';          % insert path to tiff stack here
sframe=1;						% user input: first frame to read (optional, default 1)
num2read=2000;					% user input: how many frames to read   (optional, default until the end)

Y = bigread2(nam,sframe,num2read);
if ~isa(Y,'double');    Y = double(Y);  end         % convert to double

[d1,d2,T] = size(Y);                                % dimensions of dataset
d = d1*d2;                                          % total number of pixels

%% Set parameters

K = 30;                                           % number of components to be found
tau = 4;                                          % std of gaussian kernel (size of neuron) 
p = 2;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.8;                                  % merging threshold

options = CNMFSetParms(...                      
    'd1',d1,'d2',d2,...                         % dimensions of datasets
    'search_method','ellipse','dist',3,...      % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'fudge_factor',0.98,...                      % bias correction for AR coefficients
    'merge_thr',merge_thr...                    % merging threshold
    );
%% Data pre-processing

[P,Y] = preprocess_data(Y,p);

%% fast initialization of spatial components using greedyROI and HALS

[Ain,Cin,bin,fin,center] = initialize_components(Y,K,tau,options);  % initialize

% display centers of found components
Cn =  correlation_image(Y); %max(Y,[],3); %std(Y,[],3); % image statistic (only for display purposes)
figure;imagesc(Cn);
    axis equal; axis tight; hold all;
    scatter(center(:,2),center(:,1),'mo');
    title('Center of ROIs found from initialization algorithm');
    drawnow;

%% update spatial components
Yr = reshape(Y,d,T);
clear Y;
[A,b] = update_spatial_components(Yr,Cin,fin,Ain,P,options);

%% update temporal components

[C,f,Y_res,P,S] = update_temporal_components(Yr,A,b,Cin,fin,P,options);

%% merge found components

[Am,Cm,K_m,merged_ROIs,P,Sm] = merge_components(Y_res,A,b,C,f,P,S,options);

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
        subplot(1,ln+2,ln+1); imagesc(reshape(Am(:,K_m-length(merged_ROIs)+i),d1,d2));
                title('Merged Component','fontsize',16,'fontweight','bold');axis equal; axis tight; 
        subplot(1,ln+2,ln+2);
            plot(1:T,(diag(max(C(merged_ROIs{i},:),[],2))\C(merged_ROIs{i},:))'); 
            hold all; plot(1:T,Cm(K_m-length(merged_ROIs)+i,:)/max(Cm(K_m-length(merged_ROIs)+i,:)),'--k')
            title('Temporal Components','fontsize',16,'fontweight','bold')
        drawnow;
end

%% repeat
[A2,b2] = update_spatial_components(Yr,Cm,f,Am,P,options);
[C2,f2,Y_res,P,S2] = update_temporal_components(Yr,A2,b2,Cm,f,P,options);
[C_df,~,S_df] = extract_DF_F(Yr,[A2,b2],[C2;f2],S2,K_m+1); % extract DF/F values (optional)

%% do some plotting

[A_or,C_or,S_or,P] = order_ROIs(A2,C2,S2,P);    % order components
contour_threshold = 0.95;                       % amount of energy used for each component to construct contour plot
figure;
[Coor,json_file] = plot_contours(A_or,reshape(P.sn,d1,d2),contour_threshold,1); % contour plot of spatial footprints
pause; 
%savejson('jmesh',json_file,'filename');        % optional save json file with component coordinates (requires matlab json library)
view_components(Yr,A_or,C_or,b2,f2,Cn,options);         % display all components

%% make movie

make_patch_video(A_or,C_or,b2,f2,Yr,Coor,options)
