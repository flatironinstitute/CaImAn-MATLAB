function [A,C,newcenters] = manually_refine_components(Y,A,C,centers,img,sx,options)
% function that allows to add or remove components. 
% Usage: 
% after runnning use left click to remove and right click to add components
% at specific locations. Hit enter when done
% Inputs
% ------
% centers: coordinates of the centers of the components
% img: image used to identify components (for instance NCA
% 
% Returns:
% --------
% newcenters: updated component centers

min_distance_point_selection=2;
x=1;
newcenters=[centers];
% ident_point=[0,0];
figure;
sx
while ~isempty(x)
    imagesc(img);
    axis equal; axis tight; hold all;
    scatter(newcenters(:,2),newcenters(:,1),'mo');
    title('Center of ROIs found from initialization algorithm');
%     scatter(ident_point(1),ident_point(2),'go');
    drawnow;
    [x,y,button]=ginput(1);
    disp(button)
    if ~isempty(x)
        pixel=round([x y]);
        if button==3
            disp(['Adding pixel at:' num2str(pixel)])
            newcenters=[newcenters; fliplr(pixel)];
            int_x = round(newcenters(end,1)) + (-sx:sx);
            if int_x(1)<1
                int_x = int_x + 1 - int_x(1);
            end
            if int_x(end)>options.d1
                int_x = int_x - (int_x(end)-options.d1);
            end
            int_y = round(newcenters(end,2)) + (-sx:sx);
            if int_y(1)<1
                int_y = int_y + 1 - int_y(1);
            end
            if int_y(end)>options.d2
                int_y = int_y - (int_y(end)-options.d2);
            end
            Ypatch = reshape(Y(int_x,int_y,:),(2*sx+1)^2,size(C,2));
            Ypatch = bsxfun(@minus, Ypatch, median(Ypatch,2));
            [INT_x,INT_y] = meshgrid(int_x,int_y);
            coor = sub2ind([options.d1,options.d2],INT_x(:),INT_y(:));
            Y_res = Ypatch - A(coor,:)*C;
            atemp = max(mean(Y_res,2),0);
            for i = 1:10
               ctemp = max(atemp'*Y_res,0)/norm(atemp)^2;
               atemp = max(Y_res*ctemp',0)/norm(ctemp)^2;
            end
            A(coor,end+1) = atemp/norm(atemp);
            C(end+1,:) = ctemp*norm(atemp);
            
        elseif button==1
            [m,id]=min(sum(bsxfun(@minus,newcenters,[y,x]).^2,2));            
            ident_point=[newcenters(id,2),newcenters(id,1)];
            if m<=min_distance_point_selection
                disp(['Removing point:' num2str(ident_point)]) 
                newcenters(id,:)=[]; 
                A(:,id) = [];
                C(id,:) = [];
            else
                disp('Selection too far from any point')
            end
        end
    end    
end


    function [basis, trace] = finetune2d(data, trace, nIter)
    %using matrix factorization with lasso penalty to fine-tune the basis
    %
    %Input:
    %data   M x N x T matrix, small patch containing one neuron
    %trace  initial value for trace
    %nIter  number of coordinate descent steps
    %
    %Output:
    %basis  M x N matrix, result of the fine-tuned neuron shape
    %trace  1 x T matrix, result of the neuron

        T = size(data, 3);
        if ~exist('nIter', 'var'), nIter = 1; end

        %do block coordinate descent
        for iter = 1:nIter,
            %update basis
            a = sum(trace.^2); %scale by num. of observation
            b = sum(bsxfun(@times, data, reshape(trace, [1,1,T])), 3);
            basis = max(b / a, 0);
            basisNorm = norm(basis(:));
            if basisNorm > 0, 
                basis = basis / basisNorm; 
            else
                fprintf('get degenerate basis!\n')
                break
            end

            %updating trace
            trace = bsxfun(@times, data, basis);
            trace = squeeze(sum(sum(trace, 1), 2));            
        end
    end

end