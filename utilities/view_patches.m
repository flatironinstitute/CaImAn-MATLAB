function view_patches(Y,A,C,b,f,d1,d2,sn,plot_df,make_gif)

% plot spatial components and temporal traces against filtered background
% Y:        raw data    
% A:        spatial footprints
% C:        temporal components
% b:        spatial background
% f:        temporal background
% d1:       size of FOV (rows)
% d2:       size of FOV (columns)
% sn:       normalizing constant for displaying spatial footprints (default: no normalization, other choices: P.sn)
% plot_df:  plot DF/F for temporal components (1: yes, 0: no, default: 1)
% make_gif: save a gif file with each component (1: yes, 0: no, default: 0)

% Written by:
% Eftychios A. Pnevmatikakis, Simons Foundation, 2015

T = size(C,2);
if ndims(Y) == 3
    Y = reshape(Y,d1*d2,T);
end
nr = size(A,2);     % number of ROIs
nb = size(f,1);     % number of background components
nA = full(sum(A.^2))';  % energy of each row
Y_r = spdiags(nA,0,nr,nr)\(A'*(Y- A*C - full(b)*f)) + C; 
    
if nargin < 10 || isempty(make_gif)
    make_gif = 0;
end
if nargin < 9 || isempty(plot_df)
    plot_df = 1;    
end
if nargin < 8 || isempty(sn)
    sn = ones(size(A,1),1);
end

if plot_df
    [~,Df] = extract_DF_F(Y,[A,b],[C;f],[],size(A,2)+1);
else
    Df = ones(size(A,2)+1,1);
end
    
figure;
    set(gcf,'Position',[300,300,960,480]);
for i = 1:nr+nb
    subplot(121);
    if i <= nr
        imagesc(reshape(A(:,i)./sn,d1,d2)); axis equal; axis tight;
        title(sprintf('Component %i (press any key to continue)',i),'fontsize',16,'fontweight','bold'); drawnow; %pause;
    else
        imagesc(reshape(b(:,i-nr),d1,d2)); axis equal; axis tight;
        title('Background component','fontsize',16,'fontweight','bold'); drawnow; 
    end
    subplot(122);
    if i <= nr
        plot(1:T,Y_r(i,:)/Df(i),1:T,C(i,:)/Df(i)); 
        if plot_df
            title(sprintf('Component %i (calcium DF/F value)',i),'fontsize',16,'fontweight','bold');
        else
            title(sprintf('Component %i (calcium raw value)',i),'fontsize',16,'fontweight','bold');
        end
        legend('Raw trace (filtered)','Inferred');
        xlabel('Timestep','fontsize',16,'fontweight','bold');
        drawnow; 
        if make_gif
            frame = getframe(1);
              im = frame2im(frame);
              [imind,cm] = rgb2ind(im,256);
              if i == 1;
                  imwrite(imind,cm,'results.gif','gif', 'Loopcount',inf);
              else
                  imwrite(imind,cm,'results.gif','gif','WriteMode','append');
              end
        else
            if i < nr+nb
                pause;
            end
        end
    else
        plot(1:T,f(i-nr,:)); title('Background activity','fontsize',16,'fontweight','bold');
        drawnow; 
        if make_gif
            frame = getframe(1);
              im = frame2im(frame);
              [imind,cm] = rgb2ind(im,256);
              if i == 1;
                  imwrite(imind,cm,'results.gif','gif', 'Loopcount',inf);
              else
                  imwrite(imind,cm,'results.gif','gif','WriteMode','append');
              end
        else
            if i < nr+nb
                pause;
            end
        end
    end
end