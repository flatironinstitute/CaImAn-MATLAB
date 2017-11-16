function [CC,jsf,im] = plot_contours(Aor,Cn,options,display_numbers,max_number,Coor, ln_cl, ind_show,cm)

% save and plot the contour traces of the found spatial components against
% a specified background image. The contour can be determined in two ways:
% options.thr_method = 'max': For every component pixel values below
%       options.thr_method are discarded
% options.thr_method = 'nrg': The contour is drawn around the value above
%        which a specified fraction of energy is explained (default 99%)

% INPUTS:
% Aor:              set of spatial components (matrix d x K)
% Cn:               background image (matrix d1 x d2)
% options:          options structure (optional)
% display_number:   flag for displaying the numbers of the components (optional, default: 0)
% max_number:       specify the maximum number of components to be displayed (optional, default: display all)
% Coor:             Pass contour plots to be displayed (optional)
% ln_cl:            Color of contour (default: magenta)

% OUTPUTS:
% CC:               contour plots coordinates
% jsf:              Output saved in struct format (to be saved in json format)

% Author: Eftychios A. Pnevmatikakis
%           Simons Foundation, 2016


defoptions = CNMFSetParms;
[d1,d2] = size(Cn);

if nargin < 9 || isempty(cm);
    cm=com(Aor(:,1:end),d1,d2);
end
if nargin < 8 || isempty(ind_show);
    ind_show = 1:size(Aor,2);
end

if nargin < 5 || isempty(max_number)
    max_number = size(Aor,2);
else
    max_number = min(max_number,size(Aor,2));
end
ind_show(ind_show>max_number) = [];
if nargin < 4 || isempty(display_numbers)
    display_numbers = 0;
end
if nargin < 3 || isempty(options) || isnumeric(options)
    if isnumeric(options)    % compatibility with previous version
        nrgthr = options;
        clear options
        options.thr_method = 'nrg';
        options.nrgthr = nrgthr;
        warning('plot_contours function input arguments have changed. See the file for more info.');
    else
        options = defoptions;
    end
end
if ~exist('ln_cl', 'var') || isempty(ln_cl)
    ln_cl = 'm'; % linewidth;
end
units = 'centimeters';
fontname = 'helvetica';


if ~isfield(options,'thr_method') || isempty(options.thr_method); options.thr_method = defoptions.thr_method; end
if ~isfield(options,'nrgthr') || isempty(options.nrgthr); options.nrgthr = defoptions.nrgthr; end
if ~isfield(options,'maxthr') || isempty(options.maxthr); options.maxthr = defoptions.maxthr; end

fontname = 'helvetica';
    
    if options.plot_bck_image
        im =imagesc(Cn,[min(Cn(:)),max(Cn(:))]);
    end
    axis tight; axis equal; 
    posA = get(gca,'position');
    set(gca,'position',posA);
    hold on;
    
    cmap = parula(size(Aor,2)+5); %can be 'hot'
    if ~(nargin < 6 || isempty(Coor))
        CC = Coor;
        for j = 1:length(ind_show)
            i = ind_show(j);            
            if size(Coor{i},2) > 1
                cont = medfilt1(Coor{i}')';
                cont = [cont,cont(:,2)];
                plot(cont(1,2:end),cont(2,2:end),'Color',ln_cl, 'linewidth', 2); hold on;
            end            
        end
    else
        CC = cell(size(Aor,2),1);
        CR = cell(size(Aor,2),2);
        if strcmpi(options.thr_method,'nrg')
            thr = options.nrgthr;
            for i = 1:size(Aor,2)
                A_temp = full(reshape(Aor(:,i),d1,d2));
                A_temp = medfilt2(A_temp,[3,3]);
                A_temp = A_temp(:);
                [temp,ind] = sort(A_temp(:).^2,'ascend'); 
                temp =  cumsum(temp);
                ff = find(temp > (1-thr)*temp(end),1,'first');
                if ~isempty(ff)
                    CC{i} = contour(reshape(A_temp,d1,d2),[0,0]+A_temp(ind(ff)),'LineColor',ln_cl, 'linewidth', 2);
                    fp = find(A_temp >= A_temp(ind(ff)));
                    [ii,jj] = ind2sub([d1,d2],fp);
                    CR{i,1} = [ii,jj]';
                    CR{i,2} = A_temp(fp)';
                end
                hold on;
            end       
        elseif strcmpi(options.thr_method,'max')  
            thr = options.maxthr;
            for i = 1:size(Aor,2)
                A_temp = full(reshape(Aor(:,i),d1,d2));
                A_temp = medfilt2(A_temp,[3,3]);
                A_temp(A_temp<thr*max(A_temp(:))) = 0;
                BW = bwareafilt(A_temp>0,1);                
                BW2 = bwboundaries(BW);
                if ~isempty(BW2)
                    for ii = 1:length(BW2)
                        BW2{ii} = fliplr(BW2{ii});
                        plot(BW2{ii}(:,1),BW2{ii}(:,2),'Color',ln_cl, 'linewidth', 2);
                    end
                    CC{i} = BW2{1}';
                    fp = find(BW);
                    [ii,jj] = ind2sub([d1,d2],fp);
                    CR{i,1} = [ii,jj]';
                    CR{i,2} = A_temp(fp)';
                end
                hold on;
            end
        end
    end
    if display_numbers
        lbl = strtrim(cellstr(num2str((1:size(Aor,2))')));
        for i = ind_show
            text(round(cm(i,2)),round(cm(i,1)),strtrim(cellstr(num2str(i))),'color',[1,.5 ,0],'fontsize',16,'fontname',fontname,'fontweight','bold');
        end
    end
    axis off;
    if ~(nargin < 6 || isempty(Coor))
        jsf = [];
    else
        for i = 1:size(Aor,2);
            if ~isempty(CR{i,1})
                jsf(i) = struct('id',i,...
                            'coordinates',CR{i,1}',...
                            'values',CR{i,2},...
                            'bbox',[min(CR{i,1}(1,:)),max(CR{i,1}(1,:)),min(CR{i,1}(2,:)),max(CR{i,1}(2,:))],...
                            'centroid',cm(i,:));
            end
            if i == 1
                jsf = repmat(jsf,size(Aor,2),1);
            end
        end
    end    