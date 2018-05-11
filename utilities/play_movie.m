function play_movie(movs,labels,min_mov,max_mov)
    % play list (cell) of movies with labels, press key to stop playing
    % inputs:
    % movs: cell of movies
    % labels: cell of titles for movies
    % min_mov,max_mov: min and max of movie for setting clims on movies
    
    if ~iscell(movs)
        movs = {movs};
    end
    if ~exist('labels','var')
        labels = cell(size(movs));
    end
    if ~exist('min_mov','var') || ~exist('max_mov','var')
        nml = min(1e7,numel(movs{1}));
        min_mov = quantile(movs{1}(1:nml),0.001);
        max_mov = quantile(movs{1}(1:nml),1-0.001);
    end
    
    dialogBox = uicontrol('Style', 'PushButton', 'String', 'stop','Callback', 'delete(gcbf)');

    num_movs = numel(movs);
    len_movs = size(movs{1},3);

    t = 0;
    while (ishandle(dialogBox)) && t<len_movs
        t = t+1;    
        for idx_mov = 1:num_movs
            subplot(1,num_movs,idx_mov); imagesc(movs{idx_mov}(:,:,t),[min_mov,max_mov]); 
            axis square; 
            axis off; 
            colormap('gray');
            title([labels{idx_mov}, ', frame ',num2str(t)]);
        end        
        if t == 1
            set(gcf,'Position',[100,100,1.2*length(movs)*size(movs{1},2),1.2*size(movs{1},1)]);
        end
        pause(0.001);
    end
    close()