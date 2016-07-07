function [ CC, info ] = plotROIContour( A, d1, d2, plotControl )

% This code saves and plots the contour traces of the found spatial components againsts specified background image. 
% Author: Eftychios A. Pnevmatikakis, Weijian Yang and Darcy S. Peterka 

% input parameters
% A: spatial component from CNMF 
% d1, d2: dimension of the data
% plotControl.Cn: background image
% plotControl.displayLabel: if it is 1, the ROI number will be labelled on the plot
% plotControl.ROIID: id of ROIs to be plotted
% plotControl.option: 1: The contour is drawn around the value above which a specified fraction of energy is explained (default 99%)
%                     2: The contour is drawn around the exterior boundary 
% plotControl.thr: energy fraction included in the drawn contour in option 1
% plotControl.cmap: colormap of the background image
% plotControl.lineColor: line color of the drawn contour

% output parameters
% CC: outline of the contour
% info: contour information, contains ROI id, weighted contour map, boundary box, and centroid

if isfield(plotControl,'Cn') Cn=plotControl.Cn;
else Cn=[]; end

if isfield(plotControl,'displayLabel') displayLabel=plotControl.displayLabel;
else displayLabel=1; end

if isfield(plotControl,'ROIID') ROIID=plotControl.ROIID;
else ROIID=1:size(A,2); end

if isfield(plotControl,'option') option=plotControl.option;
else option=1; end

if isfield(plotControl,'thr') thr=plotControl.thr;
else thr=0.99; end

if isfield(plotControl,'cmap') cmap=plotControl.cmap;
else cmap='gray'; end

if isfield(plotControl,'lineColor') lineColor=plotControl.lineColor;
else lineColor=[1 0 0]; end

CC = cell(length(ROIID),1);    % store contour or boundary
CR = cell(length(ROIID),2);    % store the information of the cell contour

figure;
% plot background
if ~isempty(Cn)
    imagesc(Cn,[min(Cn(:)),max(Cn(:))]);
    colormap(cmap);
end
hold on;

if option==1                   % contour is medium filtered
    for idx=1:length(ROIID)
        A_temp = full(reshape(A(:,ROIID(idx)),d1,d2));
        A_temp = medfilt2(A_temp,[3,3]);
        A_temp = A_temp(:);
        [temp,ind] = sort(A_temp(:).^2,'ascend'); 
        temp =  cumsum(temp);
        ff = find(temp > (1-thr)*temp(end),1,'first');
        if ~isempty(ff)
            CC{idx} = contour(reshape(A_temp,d1,d2),[0,0]+A_temp(ind(ff)),'LineColor',lineColor);
            CC{idx}=flipud(CC{idx});
            fp = find(A_temp >= A_temp(ind(ff)));
            [ii,jj] = ind2sub([d1,d2],fp);
            CR{idx,1} = [ii,jj]';
            CR{idx,2} = A_temp(fp)';
        end
    end
else if option==2             % contour contains the exterior boundary
    for idx=1:length(ROIID)
        A_temp = full(reshape(A(:,ROIID(idx)),d1,d2));
        [B, L]=bwboundaries(A_temp); 
        if iscell(B)
            plot(B{1}(:,2), B{1}(:,1),'linewidth',1,'color',lineColor);
            CC{idx} = B{1}';
            fp = find(L==1);
            [ii,jj] = ind2sub([d1,d2],fp);
            CR{idx,1} = [ii,jj]';
            CR{idx,2} = A_temp(fp)';
        end
    end
    end
end
   
% find the centroid
cm = com(A(:,ROIID),d1,d2);
if displayLabel
    lbl = strtrim(cellstr(num2str(ROIID')));
    text(round(cm(:,2)),round(cm(:,1)),lbl,'color',[0,0,0],'fontsize',16,'fontweight','bold');
end

% construct output 
info.id=ROIID';
info.coordinates=CR(:,1);
info.values=CR(:,2);
info.bbox=zeros(length(ROIID),4);
info.centroid=cm;
for idx = 1:length(ROIID)
    info.bbox(idx,:)=[min(CR{idx,1}(1,:)),max(CR{idx,1}(1,:)),min(CR{idx,1}(2,:)),max(CR{idx,1}(2,:))];
end

% for idx = 1:length(ROIID)
%     if ~isempty(CR{idx,1})
%         info(idx) = struct('id',ROIID(idx),...
%                         'coordinates',CR{idx,1}',...
%                         'values',CR{idx,2},...
%                         'bbox',[min(CR{idx,1}(1,:)),max(CR{idx,1}(1,:)),min(CR{idx,1}(2,:)),max(CR{idx,1}(2,:))],...
%                         'centroid',cm(idx,:));
%      end
%      if idx == 1
%          info = repmat(info,length(ROIID),1);
%      end
% end

