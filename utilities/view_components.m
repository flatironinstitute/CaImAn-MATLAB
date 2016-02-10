function view_components(Y,A,C,b,f,Cn,options)

% plot spatial components and temporal traces against filtered background
% Y:        raw data    
% A:        spatial footprints
% C:        temporal components
% b:        spatial background
% f:        temporal background
% Cn:       background image (default: mean image)
% options:  options structure

% Written by:
% Eftychios A. Pnevmatikakis, Simons Foundation, 2015

defoptions = CNMFSetParms;
if nargin < 7 || isempty(options); options = defoptions; end
if ~isfield(options,'d1') || isempty(options.d1); d1 = input('What is the total number of rows? \n'); else d1 = options.d1; end          % # of rows
if ~isfield(options,'d2') || isempty(options.d2); d2 = input('What is the total number of columns? \n'); else d2 = options.d2; end          % # of columns
% if ~isfield(options,'normalize') || isempty(options.normalize); options.normalize = ones(size(A,1),1); end
%     sn = options.normalize;
if ~isfield(options,'plot_df') || isempty(options.plot_df); options.df = defoptions.plot_df; end
    plot_df = options.plot_df;
if ~isfield(options,'make_gif') || isempty(options.make_gif); options.make_gif = defoptions.make_gif; end
    make_gif = options.make_gif; 
if ~isfield(options,'save_avi') || isempty(options.save_avi); options.save_avi = defoptions.save_avi; end
    save_avi = options.save_avi; 
if ~isfield(options,'sx') || isempty(options.sx); options.sx = defoptions.sx; end
    sx = min([options.sx,floor(d1/2),floor(d2/2)]);
if ~isfield(options,'pause_time') || isempty(options.pause_time); options.pause_time = defoptions.pause_time; end
    pause_time = options.pause_time;
if isfield(options,'name') && ~isempty(options.name); 
    name = [options.name,'_components'];
else
    name = [defoptions.name,'_components'];
end  
    
T = size(C,2);
if ndims(Y) == 3
    Y = reshape(Y,d1*d2,T);
end
if nargin < 6 || isempty(Cn);
    Cn = reshape(mean(Y,2),d1,d2);
end

nA = sqrt(sum(A.^2))';
[K,~] = size(C);
A = A/spdiags(nA,0,K,K);    % normalize spatial components to unit energy
C = spdiags(nA,0,K,K)*C;

nr = size(A,2);     % number of ROIs
nb = size(f,1);     % number of background components
%nA = full(sum(A.^2))';  % energy of each row
%Y_r = spdiags(nA,0,nr,nr)\(A'*Y- (A'*A)*C - (A'*full(b))*f) + C; 
Y_r = (A'*Y- (A'*A)*C - (A'*full(b))*f) + C;

if plot_df
    [~,Df] = extract_DF_F(Y,[A,b],[C;f],[],size(A,2)+1);
else
    Df = ones(size(A,2)+1,1);
end
    
if save_avi
    vidObj = VideoWriter([name,'.avi']);
    set(vidObj,'FrameRate',1);
    open(vidObj);
end
thr = 0.9;
fig = figure;
    set(gcf,'Position',2*[300,300,960,480]);
    set(gcf,'PaperPosition',2*[300,300,960,480]);
    int_x = zeros(nr,2*sx);
    int_y = zeros(nr,2*sx);
    cm = com(A,d1,d2);
for i = 1:nr+nb
    if i <= nr
    subplot(3,2,5);
        Atemp = reshape(A(:,i),d1,d2);
        int_x(i,:) = round(cm(i,1)) + (-(sx-1):sx);
        if int_x(i,1)<1
            int_x(i,:) = int_x(i,:) + 1 - int_x(i,1);
        end
        if int_x(i,end)>d1
            int_x(i,:) = int_x(i,:) - (int_x(i,end)-d1);
        end
        int_y(i,:) = round(cm(i,2)) + (-(sx-1):sx);
        if int_y(i,1)<1
            int_y(i,:) = int_y(i,:) + 1 - int_y(i,1);
        end
        if int_y(i,end)>d2
            int_y(i,:) = int_y(i,:) - (int_y(i,end)-d2);
        end      
        Atemp = Atemp(int_x(i,:),int_y(i,:));
        imagesc(int_x(i,:),int_y(i,:),Atemp); axis square;
    end
    subplot(3,2,[1,3]);
    if i <= nr
        imagesc(2*Cn); axis equal; axis tight; axis off; hold on; 
        A_temp = full(reshape(A(:,i),d1,d2));
        A_temp = medfilt2(A_temp,[3,3]);
        A_temp = A_temp(:)/norm(A_temp(:));
        [temp,ind] = sort(A_temp(:).^2,'ascend'); 
        temp =  cumsum(temp);
        ff = find(temp > (1-thr)*temp(end),1,'first');
        if ~isempty(ff)
            [~,ww] = contour(reshape(A_temp,d1,d2),[0,0]+A_temp(ind(ff)),'LineColor','k');
            ww.LineWidth = 2;
        end
        title(sprintf('Component %i ',i),'fontsize',16,'fontweight','bold'); drawnow; %pause;
    else
        imagesc(reshape(b(:,i-nr),d1,d2)); axis equal; axis tight;
        title('Background component','fontsize',16,'fontweight','bold'); drawnow; 
    end
    subplot(3,2,[2,4,6]);
    if i <= nr
        plot(1:T,Y_r(i,:)/Df(i),'linewidth',2); hold all; plot(1:T,C(i,:)/Df(i),'linewidth',2);
        if plot_df
            title(sprintf('Component %i (calcium DF/F value)',i),'fontsize',16,'fontweight','bold');
        else
            title(sprintf('Component %i (calcium raw value)',i),'fontsize',16,'fontweight','bold');
        end
        leg = legend('Raw trace (filtered)','Inferred');
        set(leg,'FontSize',14,'FontWeight','bold');
        drawnow; 
        hold off;
        if make_gif
            frame = getframe(fig); %getframe(1);
              im = frame2im(frame);
              [imind,clm] = rgb2ind(im,256);
              if i == 1;
                  imwrite(imind,clm,[name,'.gif'],'gif', 'Loopcount',inf);
              else
                  imwrite(imind,clm,[name,'.gif'],'gif','WriteMode','append');
              end
        else
            if i < nr+nb && ~save_avi
                fprintf('component %i. Press any key to continue.. \n', i);
                if pause_time == Inf;
                    pause;
                else
                    pause(pause_time);
                end
            end
        end
    else
        plot(1:T,f(i-nr,:)); title('Background activity','fontsize',16,'fontweight','bold');
        drawnow; 
        if make_gif
            frame = getframe(fig); %getframe(1);
              im = frame2im(frame);
              [imind,clm] = rgb2ind(im,256);
              if i == 1;
                  imwrite(imind,clm,[name,'.gif'],'gif', 'Loopcount',inf);
              else
                  imwrite(imind,clm,[name,'.gif'],'gif','WriteMode','append');
              end
        else
            if i < nr+nb && ~save_avi
                fprintf('background component %i. Press any key to continue.. \n', i-nr);
                if pause_time == Inf;
                    pause;
                else
                    pause(pause_time);
                end
            end
        end
    end
    if save_avi  
        currFrame = getframe(fig);
        writeVideo(vidObj,currFrame);    
    else
        pause(0.05);
    end
    
end

if save_avi
    close(vidObj);
end