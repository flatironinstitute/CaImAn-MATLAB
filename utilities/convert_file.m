function convert_file(input_filename,output_format,output_filename,bin_width,FOV)

% convert input file to a different format. 
% Supported input filetypes: 'mat', 'tif', 'h5', 'raw'
% Supported output filetypes: 'tif' and 'h5' (default: h5)

if isa(input_filename,'char')
    [folder_name,file_name,ext] = fileparts(input_filename);
    ext = ext(2:end);
    if strcmpi(ext,'tif') || strcmpi(ext,'tiff');
        tiffInfo = imfinfo(input_filename);
        filetype = 'tif';
        T = length(tiffInfo);
        sizY = [tiffInfo(1).Height,tiffInfo(1).Width,T];
        Y1 = imread(input_filename,'Index',1,'Info',tiffInfo);
    elseif strcmpi(ext,'mat')
        filetype = 'mem';
        Y = matfile(input_filename,'Writable',true);
        sizY = size(input_filename,'Y');
        T = sizY(end);
        if length(sizY) == 3; Y1 = Y.Y(:,:,1); elseif length(sizY) == 4; Y1 = Y.Y(:,:,:,1); end
    elseif strcmpi(ext,'hdf5') || strcmpi(ext,'h5');
        filetype = 'h5';
        fileinfo = hdf5info(input_filename);
        data_name = fileinfo.GroupHierarchy.Datasets.Name;
        sizY = fileinfo.GroupHierarchy.Datasets.Dims;
        T = sizY(end);
        Y1 = bigread2(input_filename,1);
    elseif strcmpi(ext,'raw')
        filetype = 'raw';
        fid = fopen(Y);
        if ~exist('FOV','var');
            FOV = [512,512];
            bitsize = 2;
        end
        imsize = FOV(1)*FOV(2)*bitsize;                     % Bit size of single frame
        current_seek = ftell(fid);
        fseek(fid, 0, 1);
        file_length = ftell(fid);
        fseek(fid, current_seek, -1);
        T = file_length/imsize;
        sizY = [FOV,T];
        fclose(fid);        
        Y1 = read_raw_file(input_filename,1,1);
    end    
elseif isobject(input_filename);
    filetype = 'mem';
    sizY = size(input_filename,'Y');
    T = sizY(end);
    if length(sizY) == 3; Y1 = input_filename.Y(:,:,1); elseif length(sizY) == 4; Y1 = input_filename.Y(:,:,:,1); end
    folder_name = pwd;
    file_name = 'data';
else % array loaded in memory
    filetype = 'mat';
    sizY = size(input_filename);
    T = sizY(end);
    folder_name = pwd;
    file_name = 'data';
    if length(sizY) == 3; Y1 = Y(:,:,1); elseif length(sizY) == 4; Y1 = Y(:,:,:,1); end
end

data_type = class(Y1);
if ~exist('output_format','var');
    output_format = 'h5';
end

if strcmpi(output_format,'hdf5'); output_format = 'h5'; end
if strcmpi(output_format,'tiff'); output_format = 'tif'; end

if ~exist('output_filename','var')
    output_filename = fullfile(folder_name,[file_name,'.',output_format]);
end

nd = length(sizY)-1;                          % determine whether imaging is 2d or 3d
sizY = sizY(1:nd);

if strcmpi(filetype,output_format)
    warning('\nInput file type is the same as the desired output file type. No conversion is performed');
else
    if ~exist('bin_width','var'); bin_width = round(512*512/prod(sizY)*4e3); end
    switch lower(output_format)   
        case {'hdf5','h5'}
             if exist(output_filename,'file')
                [pathstr,fname,ext] = fileparts(output_filename);             
                new_filename = fullfile(pathstr,[fname,'_',datestr(now,30),ext]);
                warning_msg = ['File ',output_filename,'already exists. Saving motion corrected file as',new_filename];            
                warning('%s',warning_msg);
                output_filename = new_filename;
            end       
            if nd == 2
                h5create(output_filename,'/mov',[sizY,Inf],'Chunksize',[sizY,bin_width],'Datatype',data_type);
            elseif nd == 3
                h5create(output_filename,'/mov',[sizY,Inf],'Chunksize',[sizY,bin_width],'Datatype',data_type);
            end
        case {'tif','tiff'}       
            opts_tiff.append = true;
            opts_tiff.big = true;
            if nd == 3
                error('Saving volumetric tiff stacks is currently not supported. Use a different filetype');
            end
        otherwise
            error('This filetype is currently not supported')
    end   
        
    for t = 1:bin_width:T
        switch filetype
            case 'tif'
                Ytm = zeros(sizY(1),sizY(2),min(t+bin_width-1,T)-t+1,data_type);
                for tt = 1:min(t+bin_width-1,T)-t+1
                    Ytm(:,:,tt) = single(imread(input_filename,'Index',t+tt-1,'Info',tiffInfo));
                end
            case 'h5'
                Ytm = h5read(input_filename,data_name,[ones(1,nd),t],[sizY(1:nd),min(t+bin_width-1,T)-t+1]);
            case 'mem'
                if nd == 2; Ytm = input_filename.Y(:,:,t:min(t+bin_width-1,T)); end
                if nd == 3; Ytm = input_filename.Y(:,:,:,t:min(t+bin_width-1,T)); end
            case 'mat'
                if nd == 2; Ytm = input_filename(:,:,t:min(t+bin_width-1,T)); end
                if nd == 3; Ytm = input_filename(:,:,:,t:min(t+bin_width-1,T)); end
            case 'raw'
                Ytm = single(read_raw_file(Y,t,min(t+bin_width-1,T)-t+1,FOV,bitsize));                
        end
        switch lower(output_format)   
            case {'hdf5','h5'}
                if nd == 2; h5write(output_filename,'/mov',Ytm,[ones(1,nd),t],[sizY(1:nd),size(Ytm,3)]); end
                if nd == 3; h5write(output_filename,'/mov',Ytm,[ones(1,nd),t],[sizY(1:nd),size(Ytm,4)]); end
            case {'tif','tiff'}
                saveastiff(cast(Ytm,data_type),output_filename,opts_tiff);
        end 
    end
end