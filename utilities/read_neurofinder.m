function imgs = read_neurofinder(foldername,output_name)
%% import neurofinder tiff multifile format

format compact
opts_tiff.append = false;
opts_tiff.big = true;
files = subdir(fullfile(foldername,'*.tif*'));
numFiles = length(files)
data_type = class(imread(files(1).name));
imgs = zeros([size(imread(files(1).name)),numFiles],data_type);
for i = 1:length(files)
    fname = files(i).name;
	imgs(:,:,i) = imread(fname);
    if mod(i,100)==0
        disp(i)
    end
end
if ~isempty(output_name)
    disp('saving to tif')
    saveastiff(cast(imgs,data_type),output_name,opts_tiff);
end
