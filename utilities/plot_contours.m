function [CC,jsf] = plot_contours(Aor,Cn,thr,display_numbers,max_number,Coor)

% save and plot the contour traces of the found spatial components againsts
% specified background image. The contour is drawn around the value above
% which a specified fraction of energy is explained (default 99%)

if nargin < 5 || isempty(max_number)
    max_number = size(Aor,2);
else
    max_number = min(max_number,size(Aor,2));
end
if nargin < 4 || isempty(display_numbers)
    display_numbers = 0;
end
if nargin < 3 || isempty(thr)
    thr = 0.995;
end

units = 'centimeters';
fontname = 'helvetica';

%fig3 = figure;
%     set(gcf, 'PaperUnits', units,'Units', units)           
%     set(gcf, 'PaperPosition',[5, 5, 12, 12])
%     set(gcf, 'Position',3*[5, 5, 12, 12])
    [d1,d2] = size(Cn);
    imagesc(Cn,[min(Cn(:)),max(Cn(:))]);
    axis tight; axis equal; 
    %set(gca,'XTick',[],'YTick',[]);
    posA = get(gca,'position');
    set(gca,'position',posA);
    %cbar = colorbar('south','TickDirection','out');
    if (0)
        cbar = colorbar('TickDirection','out');
        cpos = get(cbar,'position');
        %cpos = [posA(1),posA(2)-cpos(4)-0.01,posA(3),cpos(4)];
        ylabel(cbar,'Average neighbor correlation');
        set(cbar,'position',cpos,'TickDirection','in');
        set(cbar,'fontweight','bold','fontsize',14,'fontname',fontname);
    end
    %hold on; scatter(cm(:,2),cm(:,1),'ko'); hold off; 
    %v = axis;
    %handle = title('Correlation image and identified spatial footprints','fontweight','bold','fontsize',14,'fontname',fontname);
    hold on;
    
    cmap = hot(3*size(Aor,2));
    if ~(nargin < 6 || isempty(Coor))
        CC = Coor;
        for i = 1:size(Aor,2)
            cont = medfilt1(Coor{i}')';
            if size(cont,2) > 1
                plot(cont(1,2:end),cont(2,2:end),'Color',cmap(i+size(Aor,2),:)); hold on;
            end
        end
    else
        CC = cell(size(Aor,2),1);
        CR = cell(size(Aor,2),2);
        for i = 1:size(Aor,2)
            A_temp = full(reshape(Aor(:,i),d1,d2));
            A_temp = medfilt2(A_temp,[3,3]);
            A_temp = A_temp(:);
            [temp,ind] = sort(A_temp(:).^2,'ascend'); 
            temp =  cumsum(temp);
            ff = find(temp > (1-thr)*temp(end),1,'first');
            if ~isempty(ff)
                CC{i} = contour(reshape(A_temp,d1,d2),[0,0]+A_temp(ind(ff)),'LineColor',cmap(i+size(Aor,2),:));
                fp = find(A_temp >= A_temp(ind(ff)));
                [ii,jj] = ind2sub([d1,d2],fp);
                CR{i,1} = [ii,jj]';
                CR{i,2} = A_temp(fp)';
            end
            hold on;
        end
    end
    cm = com(Aor(:,1:end),d1,d2);
    if display_numbers
        lbl = strtrim(cellstr(num2str((1:size(Aor,2))')));
        text(round(cm(1:max_number,2)),round(cm(1:max_number,1)),lbl(1:max_number),'color',[0,0,0],'fontsize',16,'fontname',fontname,'fontweight','bold');
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