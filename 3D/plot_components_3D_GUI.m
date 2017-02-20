function plot_components_3D_GUI(Y,A,C,b,f,Cn,options)

% GUI for plotting components for the case of 3D volumetric imaging

memmaped = isobject(Y);
defoptions = CNMFSetParms;
if nargin < 7 || isempty(options); options = []; end
if ~isfield(options,'d1') || isempty(options.d1); d1 = input('What is the total number of rows? \n'); else d1 = options.d1; end          % # of rows
if ~isfield(options,'d2') || isempty(options.d2); d2 = input('What is the total number of columns? \n'); else d2 = options.d2; end          % # of columns
if ~isfield(options,'d3') || isempty(options.d2); d3 = input('What is the total number of z-stacks? \n'); else d3 = options.d3; end          % # of columns
if ~isfield(options,'plot_df') || isempty(options.plot_df); options.df = defoptions.plot_df; end
plot_df = options.plot_df;
if ~isfield(options,'full_A') || isempty(options.full_A); full_A = defoptions.full_A; else full_A = options.full_A; end

T = size(C,2);
if ~ismatrix(Y)
    Y = reshape(Y,d1*d2*d3,T);
end
if nnz(A)/numel(A) > 0.3; A = full(A); end
center = com(A,d1,d2,d3);
if size(center,2) == 2
    center(:,3) = 1;
end
center = round(center);
b = double(b);
C = double(C);
f = double(f);
nA = full(sqrt(sum(A.^2))');
[K,~] = size(C);
A = A/spdiags(nA,0,K,K);    % normalize spatial components to unit energy
C = bsxfun(@times,C,nA(:));

nb = size(f,1);     % number of background components

step = 5e3;
if memmaped
    AY = zeros(K,T);
    for j = 1:step:d
        AY = AY + A(j:min(j+step-1,d),:)'*double(Y.Yr(j:min(j+step-1,d),:));
    end
else
    if issparse(A) && isa(Y,'single')          
        if full_A
            AY = full(A)'*Y;
        else
            AY = A'*double(Y);
        end
    else
        AY = A'*Y;
    end
end
Y_r = (AY- (A'*A)*C - full(A'*double(b))*f) + C;

if plot_df
    [~,Df] = extract_DF_F(Y,[A,double(b)],[C;f],[],options);
else
    Df = ones(size(A,2)+1,1);
end

fig = figure('Visible','off');
set(gcf,'Position',2*[300,300,960,480]);
set(gcf,'PaperPosition',2*[300,300,960,480]);

thr = 0.95;
minC = max(squeeze(min(min(Cn,[],1),[],2)),0);
maxC = squeeze(max(max(Cn,[],1),[],2));

% Create a figure and axes
% ax = axes('Units','DF/F');
% Create slider

sld = uicontrol('Style', 'slider',...
    'Min',1,'Max',K+nb,'Value',1,'SliderStep',[1/(K+nb-1) 1],...
    'Position', [150 20 800 20],...
    'Callback', @surfzlim);

% Add a text uicontrol to label the slider.
txt = uicontrol('Style','text',...
    'Position',[400 45 120 20],...
    'String','Component');

% Make figure visble after adding all components
fig.Visible = 'on';
plot_component(1)

% This code uses dot notation to set properties.
% Dot notation runs in R2014b and later.
% For R2014a and earlier: set(f,'Visible','on');

    function surfzlim(source,callbackdata)
        i = source.Value;
        plot_component(round(i))
        % For R2014a and earlier:
        % i = get(source,'Value');               
    end

    function plot_component(i)
        if i <= K
            %subplot(3,2,5);
             subplot(3,3,7);
            cla
            Atemp = reshape(full(A(:,i)),d1,d2,d3);
            data = smooth3(Atemp);
            patch(isocaps(data,0.25*max(data(:))),...
               'FaceColor','interp','EdgeColor','none');
            p1 = patch(isosurface(data,0.25*max(data(:))),...
               'FaceColor','blue','EdgeColor','none');
            isonormals(data,p1)
            view(3); 
            axis vis3d tight
            camlight left; 
            %colormap jet
            lighting gouraud
            title(sprintf('Component %i',i),'fontsize',16,'fontweight','bold'); drawnow;
        end
        if i <= K
            %subplot(3,2,5);
             subplot(3,3,[2,5,8]);
            %cla
            Atemp = reshape(full(A(:,i)),d1,d2,d3);
            data = smooth3(Atemp);
            patch(isocaps(data,0.25*max(data(:))),...
               'FaceColor','interp','EdgeColor','none');
            p1 = patch(isosurface(data,0.25*max(data(:))),...
               'FaceColor','blue','EdgeColor','none');
            isonormals(data,p1)
            view(3); 
            axis vis3d %tight
            set(gca,'XLim',[1,d1],'YLim',[1,d2],'ZLim',[1,d3]);
            camlight left; 
            %colormap jet
            lighting gouraud
            title(sprintf('Components 1 - %i',i),'fontsize',16,'fontweight','bold'); drawnow;
        end        
        %subplot(3,2,[1,3]);
        subplot(3,3,[1,4]);
        if i <= K
            cla
            %,[minC(center(i,3)),maxC(center(i,3))]
            imagesc(Cn(:,:,center(i,3)),[minC(center(i,3)),maxC(center(i,3))]); axis equal; axis tight; axis off; hold on; colorbar;
            A_temp = reshape(full(A(:,i)),d1,d2,d3);
            A_temp = A_temp(:,:,center(i,3));
            A_temp = medfilt2(A_temp,[3,3]);            
            A_temp = A_temp(:);
            [temp,ind] = sort(A_temp(:).^2,'ascend');
            temp =  cumsum(temp);
            ff = find(temp > (1-thr)*temp(end),1,'first');
            if ~isempty(ff)
                [~,ww] = contour(reshape(A_temp,d1,d2),[0,0]+A_temp(ind(ff)),'LineColor','k');
                ww.LineWidth = 2;
            end
            title(sprintf('Component %i, z-plane %i',i,center(i,3)),'fontsize',16,'fontweight','bold'); drawnow; %pause;
        else
            cla
            b_temp = reshape(b(:,i-K),d1,d2,d3);
            imagesc(b_temp(:,:,ceil(d3/2))); axis equal; axis tight;
            title(['Background component, z-plane',num2str(ceil(d3/2))],'fontsize',16,'fontweight','bold'); drawnow;
        end
        %subplot(3,2,[2,4,6]);
        subplot(3,3,[3,6,9]);
        if i <= K
            plot(1:T,Y_r(i,:)/Df(i),'linewidth',2); hold all; plot(1:T,C(i,:)/Df(i),'linewidth',2);
            if plot_df
                title(sprintf('Component %i (calcium DF/F value)',i),'fontsize',16,'fontweight','bold');
            else
                title(sprintf('Component %i (calcium raw value)',i),'fontsize',16,'fontweight','bold');
            end
            leg = legend('Raw trace (filtered)','Inferred');
            set(leg,'FontSize',14,'FontWeight','bold');
            title(sprintf('Component %i',i),'fontsize',16,'fontweight','bold'); drawnow; %pause;
            drawnow;
            hold off;
        else
            plot(1:T,f(i-K,:)); title('Background activity','fontsize',16,'fontweight','bold');
            drawnow;
        end 
    end
end