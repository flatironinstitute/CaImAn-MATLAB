clear;
%% load file (courtesy of C. Lacefield and R. Bruno, Columbia University)

addpath(genpath('../../ca_source_extraction'));
nam = 'dendritic_demo.tif';
if ~exist(nam,'file')  % download file if it doesn't exist in the directory
    url = 'https://www.dropbox.com/s/wykkq4lwoo85ml5/dendritic_demo.tif.zip?dl=1';
    filename = 'dendritic_demo.zip';
    outfilename = websave(filename,url);
    unzip(filename);
end

Y = bigread2(nam);
Y = Y - min(Y(:)); 
if ~isa(Y,'double');    Y = double(Y);  end         % convert to double

[d1,d2,T] = size(Y);                                % dimensions of dataset
d = d1*d2;                                          % total number of pixels

%% view data

nY = quantile(Y(:),0.01);
mY = quantile(Y(:),0.995);

figure;
for t = 1:4:T
    imagesc(Y(:,:,t),[nY,mY]);
    title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14);
    drawnow;
    pause(0.01);
end

%% Set parameters

K = 50;                                           % number of components to be found
tau = [];                                         % std of gaussian kernel (size of neuron - not needed for dendritic data) 
p = 0;                                            % No AR dynamics for dendritic data
merge_thr = 0.8;                                  % merging threshold

options = CNMFSetParms(...                      
    'd1',d1,'d2',d2,...                           % dimensions of datasets
    'init_method','HALS',...                      % initialize algorithm with plain NMF  
    'max_iter_hals_in',50,...                     % maximum number of iterations
    'search_method','dilate',...                  % search locations when updating spatial components
    'temporal_iter',2,...                         % number of block-coordinate descent steps 
    'merge_thr',0.8,...                           % merging threshold
    'conn_comp',false,...                         % do not limit to largest connected component for each found component
    'maxthr',0.05...                              % for every component set pixels that are below max_thr*max_value to 0 
    );
%% Data pre-processing

[P,Y] = preprocess_data(Y,p);

%% fast initialization of spatial components using greedyROI and HALS

[Ain,Cin,bin,fin,center] = initialize_components(Y,K,tau,options,P);  % initialize

Yr = reshape(Y,d,T);
plot_dend_components_GUI(Yr,Ain,Cin,bin,fin,options);  % view the components

    
%% update spatial components

[A,b,Cin] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);

%% update temporal components

[C,f,P,S] = update_temporal_components(Yr,A,b,Cin,fin,P,options);

%% merge found components
[Am,Cm,K_m,merged_ROIs,P,Sm] = merge_components(Yr,A,b,C,f,P,S,options);

display_merging = 1; % flag for displaying merging example
if and(display_merging, ~isempty(merged_ROIs))
    i = 1; %randi(length(merged_ROIs));
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
else
    fprintf('No components were merged. \n')
end

%% order and extract DF/F

[A_or,C_or,S_or,P] = order_ROIs(A,C,S,P); % order components
K_m = size(C_or,1);
[C_df,~] = extract_DF_F(Yr,[A_or,b],[C_or;f],K_m+1); % extract DF/F values (optional)

%% display components

plot_dend_components_GUI(Yr,A_or,C_or,b,f,options)

%% make movie

make_dendritic_video(A_or,C_or,b,f,Yr,d1,d2)
