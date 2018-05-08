function imData=read_file(path_to_file,sframe,num2read,options,im_info)

% Reads uncompressed multipage .tiff, .hdf5, .avi or .raw files 
% Usage:  my_data=read_file('path_to_data_file, start frame, num to read);

% INPUTS:
% path_to_file:     location of file to be read
% sframe:           first frame to read (optional, default: 1)
% num2read:         number of frames to read (optional, default: read the whole file)
% options:          options for reading .raw or .bin files
% im_info:          information about the file (if already present)


% OUTPUT:
% imData:           data in array format 

% Written by Eftychios A. Pnevmatikakis, Simons Foundation

if nargin<2 || isempty(sframe); sframe = 1; end
if nargin<3 || isempty(num2read); num2read = Inf; end

[~,~,ext] = fileparts(path_to_file);

if strcmpi(ext,'.tiff') || strcmpi(ext,'.tif') || strcmpi(ext,'.btf')  
    if ~exist('im_info','var')
        im_info = imfinfo(path_to_file);
    end
    TifLink = Tiff(path_to_file, 'r');
    num2read = min(num2read,length(im_info)-sframe+1);
    imData = zeros(im_info(1).Height,im_info(1).Width,num2read,'like',TifLink.read());
    for i=1:num2read
       TifLink.setDirectory(i+sframe-1);
       imData(:,:,i)=TifLink.read();
    end
    TifLink.close()
    %imData = loadtiff(path_to_file,sframe,num2read);    
elseif strcmpi(ext,'.hdf5') || strcmpi(ext,'.h5')
%     info = hdf5info(path_to_file);
%     dims = info.GroupHierarchy.Datasets.Dims;
%     name = info.GroupHierarchy.Datasets.Name;
    info = h5info(path_to_file);
    dims = info.Datasets.Dataspace.Size;
    name = info.Datasets.Name;    
    if nargin < 3
        num2read = dims(end)-sframe+1;
    end
    num2read = min(num2read,dims(end)-sframe+1);
%    imData = h5read(path_to_file,name,[ones(1,length(dims)-1),sframe],[dims(1:end-1),num2read]);
    imData = h5read(path_to_file,['/',name],[ones(1,length(dims)-1),sframe],[dims(1:end-1),num2read]);
elseif strcmpi(ext,'.avi')
    v = VideoReader(path_to_file);
    if v.Duration*v.FrameRate < sframe
        imData = [];
    else
        if nargin < 3
            num2read = v.Duration*v.FrameRate-sframe+1;
        end
        num2read = min(num2read,v.Duration*v.FrameRate-sframe+1);
        Y1 = readFrame(v);
        imData = zeros(v.Height,v.Width,num2read,'like',Y1);
        i = 1;
        if sframe == 1
            imData(:,:,i-sframe+1) = Y1;
        end
        while hasFrame(v)  && (i - sframe + 1 < num2read)
            video = readFrame(v);
            i = i + 1;
            if i >= sframe
                imData(:,:,i-sframe+1) = video;
            end
            if i - sframe + 1 >= num2read
                break;
            end
        end
    end
elseif strcmpi(ext,'.raw')
    if nargin < 4 || ~isfield(options,'d1') || isempty(options.d1);
        options.d1 = input('What is the total number of rows? \n');
    end
    if nargin < 4 || ~isfield(options,'d2') || isempty(options.d2);
        options.d2 = input('What is the total number of columns? \n');
    end
    if nargin < 4 || ~isfield(options,'bitsize') || isempty(options.bitsize);
        options.bitsize = 2;
    end   
    imData = read_raw_file(path_to_file,sframe,num2read,[options.d1,options.d2],options.bitsize);
else
    error('Unknown file extension. Only .tiff, .avi and .hdf5 files are currently supported');
end