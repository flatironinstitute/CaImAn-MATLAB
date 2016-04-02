

% Plot spatial component and S,C
unit=zeros(d1,d2);

% Open decent size figure
figure
scrsz = get(groot,'ScreenSize');
set(gcf, 'Position',[50 50 scrsz(3)-100 scrsz(4)*2/3-150])
pan2 = panel();
pan2.pack('h', {1/4 []} )

for i=1:size(S_df,1)
    % Spatial component
    pan2(1).select()
    unit=unit + reshape(A(:,i), d1, d2);
    imagesc(unit);
    axis image; 
    
    % S and C
    pan2(2).select()
    plot([S_df(i,:); C_df(i,:)]')
    pause
end