function make_patch_video(A,C,b,f,Y,cont,options)

defoptions = CNMFSetParms;
if nargin < 7
    options = [];
end

if isfield(options,'ind');  % indeces of components to be shown
    ind = options.ind;
else
    ind = defoptions.ind;       
end

if isfield(options,'skip_frame'); % skip frames when showing the video
    skp = options.skip_frame;
else
    skp = defoptions.skip_frame;
end

if isfield(options,'sx');   % size of patch in pixels
    sx = options.sx;
else
    sx = defoptions.sx;        
end

if isfield(options,'make_avi');  % save an avi video
    make_avi = options.make_avi;
else
    make_avi = defoptions.make_avi;       
end

if isfield(options,'name');  % name of video
    name = options.name;
else
    name = defoptions.name;  
end

if isfield(options,'show_background');
    sb = options.show_background;
else
    sb = defoptions.show_background;
end

if isfield(options,'show_contours');
    sc = options.show_contours;
else
    sc = defoptions.show_contours;
end

if isfield(options,'cmap');
    cmap = options.cmap;
else
    cmap = defoptions.cmap;
end

d1 = options.d1;
d2 = options.d2;

if isempty(cont) || nargin < 6; sc = 0; end

T = size(C,2);

if ismatrix(Y)
    Y = reshape(Y,d1,d2,T);
end
sx = min(sx,floor(min(d1,d2)/2));

cm = com(A,d1,d2);

C_rec = reshape(A*C+b*f,d1,d2,T);
C_np = reshape(b*f,d1,d2,T);
if make_avi
    vidObj = VideoWriter(name);
    set(vidObj,'FrameRate',15);
    open(vidObj);
end
fig = figure; %colormap('bone');
    set(gcf, 'PaperUnits', 'inches','Units', 'inches')           
    set(gcf, 'PaperPositionMode', 'manual')
    set(gcf, 'PaperPosition',1.5*[0,0, 14, 10.5]/1.5)
    set(gcf, 'Position',1.5*[2,2, 14, 10.5]/1.5)%Yr = reshape(Y,d1,d2,T);
    colormap(cmap);
int_x = zeros(length(ind),2*sx);
int_y = zeros(length(ind),2*sx);
A_com = zeros(2*sx,2*sx,length(ind));
for i = 1:length(ind)
    int_x(i,:) = round(cm(ind(i),1)) + (-(sx-1):sx);
    if int_x(i,1)<1
        int_x(i,:) = int_x(i,:) + 1 - int_x(i,1);
    end
    if int_x(i,end)>d1
        int_x(i,:) = int_x(i,:) - (int_x(i,end)-d1);
    end
    int_y(i,:) = round(cm(ind(i),2)) + (-(sx-1):sx);
    if int_y(i,1)<1
        int_y(i,:) = int_y(i,:) + 1 - int_y(i,1);
    end
    if int_y(i,end)>d2
        int_y(i,:) = int_y(i,:) - (int_y(i,end)-d2);
    end
    A_temp = reshape(A(:,ind(i)),d1,d2);
    A_com(:,:,i) = A_temp(int_x(i,:),int_y(i,:));
end
mY = 0.7*max(Y(:));
mC = 0.7*max(C_rec(:));
mB = 0.95*max(C_np(:));

up = zeros(1,4);
for i = 1:length(ind)
    up(i) = 0.95*max(max(A_com(:,:,i)))*max(C(ind(i),:));
    ff = find(cont{ind(i)}(1,2:end)<0.1);
    cont{ind(i)}(:,ff+1) = [];
end

for t = 1:skp:T
    subplot(4,6,[1,2,7,8]); imagesc(Y(:,:,t),[0,mY]); 
        title('Raw data','fontweight','bold','fontsize',16);
        axis equal; axis tight; axis off;
        if sc
            hold on;
            for i = 1:length(ind)
                plot(cont{ind(i)}(1,2:end),cont{ind(i)}(2,2:end),'w'); hold on;
                text(min(cont{ind(i)}(1,2:end))-2,min(cont{ind(i)}(2,2:end))-2,num2str(ind(i)),'color','w','fontsize',14,'fontweight','bold'); hold on;
            end
            hold off;
        end
    subplot(4,6,[3,4,9,10]); imagesc(C_rec(:,:,t) - (1-sb)*C_np(:,:,t),[0,mC]); %title(sprintf('%i',t));
        title('Denoised','fontweight','bold','fontsize',16);
        xlabel(sprintf('Timestep %i out of %i',t,T),'fontweight','bold','fontsize',16);
        axis equal; axis tight; set(gca,'XTick',[],'YTick',[]);
        if sc
            hold on;
            for i = 1:length(ind)
                plot(cont{ind(i)}(1,2:end),cont{ind(i)}(2,2:end),'w'); hold on;
            end
            hold off;            
        end
    subplot(4,6,[5,6,11,12]); imagesc(Y(:,:,t)-C_rec(:,:,t),[-mC/2,mC/2]);
        title('Residual 2x','fontweight','bold','fontsize',16);
        axis equal; axis tight; axis off;
    subplot(4,6,[13,14,19,20]); imagesc(C_np(:,:,t),[0,mB]);
        title({'Background'; 'synchronized activity'},'fontweight','bold','fontsize',16);
        axis equal; axis tight; axis off;

    for i = 1:4        
        subplot(4,6,14+i); imagesc(A_com(:,:,i)*C(ind(i),t),[0,up(i)/1.2]); axis square; 
        if sc
            hold on;
            plot(cont{ind(i)}(1,2:end)-int_y(i,1)+1  ,cont{ind(i)}(2,2:end)-int_x(i,1)+1,'w');
            hold off;            
        end
        xlabel(sprintf('Component %i',ind(i)),'fontweight','bold','fontsize',14);
        set(gca,'XTick',[],'YTick',[]);
        if i == 2
            handle = title('Representative spatial components','fontweight','bold','fontsize',16);
            set(handle,'position',get(handle,'position')./[0.43,1,1])
        end
        subplot(4,6,20+i); imagesc(Y(int_x(i,:),int_y(i,:),t),[0,up(i)+mB/2]); axis square; axis off;
        if sc
            hold on;
            plot(cont{ind(i)}(1,2:end)-int_y(i,1)+1  ,cont{ind(i)}(2,2:end)-int_x(i,1)+1,'w');
            hold off;            
        end
        if i == 2
            handle = title('Raw data patches','fontweight','bold','fontsize',16);
            set(handle,'position',get(handle,'position')./[0.43,1,1])
        end        
    end
    drawnow;
    if make_avi  
        currFrame = getframe(fig);
        writeVideo(vidObj,currFrame);    
    else
        pause(0.05);
    end
end
if make_avi
    close(vidObj);
end

end