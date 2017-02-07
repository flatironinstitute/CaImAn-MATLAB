function writeTiff(I, filename)
%% write tiff stack 
imwrite(I(:, :, 1), filename); 
for m=2:size(I, 3)
    imwrite(I(:, :, m),  filename, 'writemode', 'append'); 
end