function newcenters=manually_refine_components(centers,img)
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
        elseif button==1
            [m,id]=min(sum(bsxfun(@minus,newcenters,[y,x]).^2,2));            
            ident_point=[newcenters(id,2),newcenters(id,1)];
            if m<=min_distance_point_selection
                disp(['Removing point:' num2str(ident_point)]) 
                newcenters(id,:)=[]; 
            else
                disp('Selection too far from any point')
            end
        end
    end    
end