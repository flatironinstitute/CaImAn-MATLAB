function play_movie(movs,labels,min_mov,max_mov)
    % play list (cell) of movies with labels, press key to stop playing
    % inputs:
    % movs: cell of movies
    % labels: cell of titles for movies
    % min_mov,max_mov: min and max of movie for setting clims on movies
    
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
        pause(0.001);
    end
    close()