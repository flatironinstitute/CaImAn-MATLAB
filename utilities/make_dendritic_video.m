function make_dendritic_video(A,C,b,f,Y,d1,d2,param)

% use that to make movies for dendritic/axonal imaging data

if nargin < 8
    param = [];
end

if isfield(param,'ind');  % indeces of components to be shown
    ind = param.ind;
else
    ind = 1:7;       
end

if isfield(param,'skip_frame'); % skip frames when showing the video
    skp = param.skip_frame;
else
    skp = 1;
end

if isfield(param,'make_avi');  % save an avi video
    make_avi = param.make_avi;
else
    make_avi = 0;       
end

if isfield(param,'name');  % # of repetitions
    name = param.name;
else
    name = ['video_',datestr(now,30),'.avi'];        % default number of repetitions
end

T = size(C,2);

if ismatrix(Y)
    Yr = Y;
    Y = reshape(Y,d1,d2,T);
else
    Yr = reshape(Y,d1*d2,T);
end

nr = size(A,2);
bk_ind = nr + (1:size(b,2));
A = [A,b];
C = [C;f];

ind = [ind,size(A,2)];
Ym = zeros(2*d1,4*d2,T);
A = A';
C = C';
cnt = 0;
for i = 1:4
    for j = 1:2
        cnt = cnt + 1;
        id = ind(cnt);
        Ar = A(id,:)/max(A(id,:));
        Cr = C(:,id)/max(C(:,id));
        Ym((j-1)*d1+(1:d1),(i-1)*d2+(1:d2),:) = reshape(Ar'*Cr',d1,d2,T);
    end
end
Ym(d1:d1:end,:,:) = 1;
Ym(:,d2:d2:end,:) = 1;
Ym(d1+1:d1:end,:,:) = 1;
Ym(:,d2+1:d2:end,:) = 1;

mmY = max(Y(:));
fig = figure; 
rm_pix = find(max(C(:,ind(1:end-1)),[],2)>0.1);
set(gcf, 'Units', 'pixels')           
set(gcf, 'Position',[100,100, 960, 960]/1.2)
if make_avi
    vidObj = VideoWriter(name);
    set(vidObj,'FrameRate',10);
    open(vidObj);
end
indeces = union(1:35,bk_ind);
Y_noise = Yr - A(indeces,:)'*C(:,indeces)';
Y_noise = reshape(Y_noise,d1,d2,T);
ind_rec = setdiff(1:35,bk_ind);
Y_rec =  A(ind_rec,:)'*C(:,ind_rec)';
Y_rec = reshape(Y_rec,d1,d2,T);
cbm = max(Y_rec(:));
nR = quantile(Y_rec(:),0.02);
mR = quantile(Y_rec(:),0.999);
mY = quantile(Y(:),0.99);
nN = quantile(Y_noise(:),0.02);
mN = quantile(Y_noise(:),0.99);
for t = 1:skp:length(rm_pix)
    i = rm_pix(t);
    subplot(3,3,1); 
        imagesc(Y(:,:,i),[nR,mY]); axis equal; axis tight; 
        posf1 = get(gca,'position');
        colorbar('location','southoutside');
        set(gca,'position',posf1);
        title('Raw Data','Fontsize',16,'Fontweight','bold');
        set(gca,'XTick',[],'YTick',[]);
        pos = get(gca,'Position') + [-0.02,-0.025,0.05,0.05];
        set(gca,'Position',pos);

    subplot(3,3,2);
        imagesc(Y_rec(:,:,i),[nR,mR]); axis equal; axis tight; 
        title('Denoised Data','Fontsize',16,'Fontweight','bold');
        set(gca,'XTick',[],'YTick',[]);
        set(gca,'XTick',[],'YTick',[]);
        xlabel(sprintf('Timestep %i',i),'Fontsize',16,'Fontweight','bold');
        pos = get(gca,'Position') + [-0.01,-0.025,0.05,0.05];
        set(gca,'Position',pos);
        
    subplot(3,3,3); 
        imagesc(Y_noise(2:end-1,2:end-1,i),[nN,mN]); axis equal; axis tight;
        title('Residual','Fontsize',16,'Fontweight','bold'); axis off;  
        pos = get(gca,'Position') + [-0.01,-0.025,0.05,0.05];
        hc = colorbar('location','southoutside');
        set(gca,'position',pos);        
        
    subplot(3,3,4:9); imagesc(Ym(2:end-1,2:end-1,i),[0,1]); axis equal; axis tight; axis off;
    if t == 1
        pos_full = get(gca,'Position');
        pos_full = pos_full + [-0.065,-0.05,0.12,0.14];
    end
    title('Extracted Components','Fontsize',16,'Fontweight','bold');
    set(gca,'XTick',[],'YTick',[],'Position',pos_full);
    drawnow; pause(0.01);
    if make_avi
        currFrame = getframe(fig);
        writeVideo(vidObj,currFrame);
    end
    if t<T; 
        clf; 
    else  
        if make_avi
            currFrame = getframe(fig);
            writeVideo(vidObj,currFrame);
        end
    end    
end

if make_avi
    close(vidObj);
end