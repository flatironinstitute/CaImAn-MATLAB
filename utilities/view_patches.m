function view_patches(Y,A,C,b,f,d1,d2,sn,ind_neur)

T = size(C,2);
if ndims(Y) == 3
    Y = reshape(Y,d1*d2,T);
end
nr = size(A,2);     % number of ROIs
nb = size(f,1);     % number of background components
nA = full(sum(A.^2))';  % energy of each row
Y_r = spdiags(nA,0,nr,nr)\(A'*(Y-full(b)*f)); 
    

if nargin < 9
    ind_neur = ones(nr,1);
end
if nargin < 8
    sn = ones(size(A,1),1);
end

make_gif = 0;

figure;
    set(gcf,'Position',[300,300,960,480]);
for i = 1:nr+nb
    subplot(121);
    if i <= nr
        imagesc(reshape(A(:,i)./sn,d1,d2)); axis equal; axis tight;
        title({sprintf('Plane %i',ind_neur(i));sprintf('Component %i (press any key to continue)',i)},'fontsize',16,'fontweight','bold'); drawnow; %pause;
    else
        imagesc(reshape(b(:,i-nr),d1,d2)); axis equal; axis tight;
        title('Background component','fontsize',16,'fontweight','bold'); drawnow; 
    end
    subplot(122);
    if i <= nr
        plot(1:T,Y_r(i,:),1:T,C(i,:)); 
        title(sprintf('Component %i (calcium)',i),'fontsize',16,'fontweight','bold');
        legend('Raw trace (filtered)','Inferred');
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