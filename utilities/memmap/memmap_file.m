function data = memmap_file(filename,sframe,num2read,chunksize)

% read a stacked tiff array, reshapes it to 2d array and saves it a mat file that can be memory mapped. 
% The file is saved both in its original format and as a reshaped 2d matrix

% INPUTS
% filename:     path to tiff file
% sframe:       (optional) first frame to read (default: 1)
% num2read:     (optional) number of frames to read (default: until the end of file)
% chunksize:    (optional) read the file in chunks for big files (default: read all at once)

% OUTPUT
% data:         object with the data containing:
%   Yr:         reshaped data file
% sizY:         dimensions of original size
%   nY:         minimum value of Y (not subtracted from the dataset)

% Author: Eftychios A. Pnevmatikakis, Simons Foundation, 2016

if nargin < 4 
    chunksize = [];
end

if nargin < 3
    num2read = [];
end

if nargin < 2 || isempty(sframe)
    sframe = 1;
end

if isempty(chunksize)
    Y = bigread2(filename,sframe,num2read);
    sizY = size(Y);
    Yr = reshape(Y,prod(sizY(1:end-1)),[]);
    nY = min(Yr(:));
    %Yr = Yr - nY;
    save([filename(1:end-3),'mat'],'Yr','Y','nY','sizY','-v7.3');
    data = matfile([filename(1:end-3),'mat'],'Writable',true);
else
    info = imfinfo(filename);
    blah=size(info);
    if isempty(num2read)
        numFrames = blah(1);
    else
        numFrames = min(blah(1),sframe+num2read-1);
    end
    nY = Inf;
    data = matfile([filename(1:end-3),'mat'],'Writable',true);
    for i = sframe:chunksize:numFrames
        Ytemp = bigread2(filename,i,min(chunksize,numFrames-i+1));
        sizY = size(Ytemp);
        Yr = reshape(Ytemp,prod(sizY(1:end-1)),[]);
        nY = min(nY,min(Yr(:)));
        data.Yr(1:prod(sizY(1:end-1)),i-1+(1:size(Yr,2))) = Yr;
        if length(sizY) == 3
            data.Y(1:sizY(1),1:sizY(2),i-1+(1:size(Yr,2))) = Ytemp;
        elseif legnth(sizY) == 4
            data.Y(1:sizY(1),1:sizY(2),1:sizY(3),i-1+(1:size(Yr,2))) = Ytemp;
        end
    end
    data.nY = nY;
    data.sizY = [sizY(1:end-1),numFrames-sframe+1];
end