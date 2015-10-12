function make_patch_video(A,C,b,f,Y,d1,d2,param)

if nargin < 8
    param = [];
end

    if isfield(param,'ind');  % indeces of components to be shown
        ind = param.ind;
    else
        ind = 1:4;       
    end

    if isfield(param,'skip_frame'); % skip frames when showing the video
        skp = param.skip_frame;
    else
        skp = 1;
    end
    
    if isfield(param,'sx');   % size of patch in pixels
        sx = param.sx;
    else
        sx = 32;        
    end

    if isfield(param,'make_avi');  % save an avi video
        make_avi = param.make_avi;
    else
        make_avi = 0;       
    end

    if isfield(param,'name');  % name of video
        name = param.name;
    else
        name = ['video_',datestr(now,30),'.avi'];  
    end

    if isfield(param,'show_background');
        sb = param.show_background;
    else
        sb = 1;
    end
    
    T = size(C,2);

    if ismatrix(Y)
        Y = reshape(Y,d1,d2,T);
    end

    Coor.x = kron(ones(d2,1),(1:d1)');
    Coor.y = kron((1:d2)',ones(d1,1));
    nr = size(A,2);
    cm = zeros(nr,2);  % vector for center of mass							   
    cm(:,1) = Coor.x'*A./sum(A);
    cm(:,2) = Coor.y'*A./sum(A);

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
        colormap('bone');
    int_x = zeros(4,2*sx);
    int_y = zeros(4,2*sx);
    A_com = zeros(2*sx,2*sx,4);
    for i = 1:4
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
    for i = 1:4
        up(i) = 0.95*max(max(A_com(:,:,i)))*max(C(ind(i),:));
    end

    for t = 1:skp:T
        subplot(4,6,[1,2,7,8]); imagesc(Y(:,:,t),[0,mY]); 
        title('Raw data','fontweight','bold','fontsize',16);
        axis square; axis off;
        subplot(4,6,[3,4,9,10]); imagesc(C_rec(:,:,t) - (1-sb)*C_np(:,:,t),[0,mC]); %title(sprintf('%i',t));
            title('Denoised','fontweight','bold','fontsize',16);
            xlabel(sprintf('Timestep %i out of %i',t,T),'fontweight','bold','fontsize',16);
            axis square; set(gca,'XTick',[],'YTick',[]);
        subplot(4,6,[5,6,11,12]); imagesc(Y(:,:,t)-C_rec(:,:,t),[-mC/2,mC/2]);
            title('Residual 2x','fontweight','bold','fontsize',16);
            axis square; axis off;
        subplot(4,6,[13,14,19,20]); imagesc(C_np(:,:,t),[0,mB]);
            title({'Background'; 'synchronized activity'},'fontweight','bold','fontsize',16);
            axis square; axis off;

        for i = 1:4        
            subplot(4,6,14+i); imagesc(A_com(:,:,i)*C(ind(i),t),[0,up(i)/1.2]); axis square; 
            xlabel(sprintf('Component %i',ind(i)),'fontweight','bold','fontsize',14);
            set(gca,'XTick',[],'YTick',[]);
            if i == 2
                handle = title('Representative spatial components','fontweight','bold','fontsize',16);
                set(handle,'position',get(handle,'position')./[0.43,1,1])
            end
            subplot(4,6,20+i); imagesc(Y(int_x(i,:),int_y(i,:),t),[0,up(i)+mB/2]); axis square; axis off;
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