function imData=read_file_old(path_to_file,sframe,num2read)

% Reads uncompressed multipage .tiff, .hdf5, or .avi files 
% Usage:  my_data=read_file('path_to_data_file, start frame, num to read);

% INPUTS:
% path_to_file:     location of file to be read
% sframe:           first frame to read (optional, default: 1)
% num2read:         number of frames to read (optional, default: read the whole file)

% OUTPUT:
% imData:           data in array format 

% Written by Eftychios A. Pnevmatikakis, Simons Foundation

if nargin<2 || isempty(sframe); sframe = 1; end

[~,~,ext] = fileparts(path_to_file);

if strcmpi(ext,'.tiff') || strcmpi(ext,'.tif');
    
    %get image info
    tiffInfo = imfinfo(path_to_file);
    T = length(tiffInfo);
    if nargin < 3 || isempty(num2read); num2read = T - sframe + 1; end
    num2read = min(num2read,T - sframe + 1);
    
    Y1 = imread(path_to_file,'Index',sframe,'Info',tiffInfo);
    imData = zeros([size(Y1),num2read],'like',Y1);
    nd = ndims(Y1);
    if nd == 2
        imData(:,:,1) = Y1;   
        for t = sframe+1:sframe+num2read-1
            imData(:,:,t-sframe+1) = imread(path_to_file,'Index',t,'Info',tiffInfo);
        end
    elseif nd == 3
        imData(:,:,:,1) = Y1;   
        for t = sframe+1:sframe+num2read-1
            imData(:,:,:,t-sframe+1) = imread(path_to_file,'Index',t,'Info',tiffInfo);
        end        
    end
    
elseif strcmpi(ext,'.hdf5') || strcmpi(ext,'.h5');
    info = hdf5info(path_to_file);
    dims = info.GroupHierarchy.Datasets.Dims;
    name = info.GroupHierarchy.Datasets.Name;
    if nargin < 3
        num2read = dims(end)-sframe+1;
    end
    num2read = min(num2read,dims(end)-sframe+1);
    imData = h5read(path_to_file,name,[ones(1,length(dims)-1),sframe],[dims(1:end-1),num2read]);
elseif strcmpi(ext,'.avi')
    v = VideoReader(path_to_file);
    if nargin < 3
        num2read = v.Duration*v.FrameRate-sframe+1;
    end
    Y1 = readFrame(v);
    imData = zeros(v.Height,v.Width,num2read,'like',Y1);
    i = 0;
    while hasFrame(v)
        video = readFrame(v);
        i = i + 1;
        if i >= sframe
            imData(:,:,i-sframe+1) = video;
        end
        if i - sframe + 1 >= num2read
            break;
        end
    end
else
    error('Unknown file extension. Only .tiff, .avi and .hdf5 files are currently supported');
end