function imData=bigread2(path_to_file,sframe,num2read)
%reads tiff files in Matlab bigger than 4GB, allows reading from sframe to sframe+num2read-1 frames of the tiff - in other words, you can read page 200-300 without rading in from page 1.
%based on a partial solution posted on Matlab Central (http://www.mathworks.com/matlabcentral/answers/108021-matlab-only-opens-first-frame-of-multi-page-tiff-stack)
%Darcy Peterka 2014, v1.0
%Darcy Peterka 2014, v1.1 
%Darcy Peterka 2016, v1.2(bugs to dp2403@columbia.edu)
%Eftychios Pnevmatikakis 2016, v1.3 (added hdf5 support)
%Program checks for bit depth, whether int or float, and byte order.  Assumes uncompressed, non-negative (i.e. unsigned) data.
%
% Usage:  my_data=bigread('path_to_data_file, start frame, num to read);
% "my_data" will be your [M,N,frames] array.
%Will do mild error checking on the inputs - last two inputs are optional -
%if nargin == 2, then assumes second is number of frames to read, and
%starts at frame 1

[~,~,ext] = fileparts(path_to_file);

if strcmpi(ext,'.tiff') || strcmpi(ext,'.tif');
    
    %get image info
    info = imfinfo(path_to_file);

    % if ~isfield(info,'ImageDescription')
    %    blah=size(info);
    %    numFrames= blah(1);
    % else
    %     he=info.ImageDescription;
    %     numFramesStr = regexp(he, 'images=(\d*)', 'tokens');
    %     numFrames = str2double(numFramesStr{1}{1});
    % end
    blah=size(info);
    numFrames= blah(1);

    num_tot_frames=numFrames;

    %should add more error checking for args... very ugly code below.  works
    %for me after midnight though...
    if nargin<2
        sframe = 1;
    end
    if nargin<3
        num2read=numFrames-sframe+1;
    end
    if sframe<=0
        sframe=1;
    end
    if num2read<1
        num2read=1;
    end
    if sframe>num_tot_frames
        sframe=num_tot_frames;
        num2read=1;
        display('starting frame has to be less than number of total frames...');
    end
    if (num2read+sframe<= num_tot_frames+1)
        lastframe=num2read+sframe-1;
    else
        num2read=numFrames-sframe+1;
        lastframe=num_tot_frames;
        display('Hmmm...just reading from starting frame until the end');
    end


    bd=info.BitDepth;
    he=info.ByteOrder;
    bo=strcmp(he,'big-endian');
    if (bd==64)
        form='double';
    elseif(bd==32)
        form='single'
    elseif (bd==16)
        form='uint16';
    elseif (bd==8)
        form='uint8';
    end


    % Use low-level File I/O to read the file
    fp = fopen(path_to_file , 'rb');
    % The StripOffsets field provides the offset to the first strip. Based on
    % the INFO for this file, each image consists of 1 strip.

    he=info.StripOffsets;
    %finds the offset of each strip in the movie.  Image does not have to have
    %uniform strips, but needs uniform bytes per strip/row.
    idss=max(size(info(1).StripOffsets));
    ofds=zeros(numFrames);
    for i=1:numFrames
        ofds(i)=info(i).StripOffsets(1);
        %ofds(i)
    end
    sframemsg = ['Reading from frame ',num2str(sframe),' to frame ',num2str(num2read+sframe-1),' of ',num2str(num_tot_frames), ' total frames'];
    disp(sframemsg)
    pause(.2)
    %go to start of first strip
    fseek(fp, ofds(1), 'bof');
    %framenum=numFrames;
    framenum=num2read;
    imData=cell(1,framenum);

    he_w=info.Width;
    he_h=info.Height;
    % mul is set to > 1 for debugging only
    mul=1;
    if strcmpi(form,'uint16') || strcmpi(form,'uint8')
        if(bo)
            for cnt = sframe:lastframe
                %cnt;
                fseek(fp,ofds(cnt),'bof');
                tmp1 = fread(fp, [he_w he_h*mul], form, 0, 'ieee-be')';
                imData{cnt-sframe+1}=cast(tmp1,form);
            end
        else
            for cnt=sframe:lastframe
                % cnt;
                fseek(fp,ofds(cnt),'bof');
                tmp1 = fread(fp, [he_w he_h*mul], form, 0, 'ieee-le')';
                imData{cnt-sframe+1}=cast(tmp1,form);
            end
        end
    elseif strcmpi(form,'single')
        if(bo)
            for cnt = sframe:lastframe
                %cnt;
                fseek(fp,ofds(cnt),'bof');
                tmp1 = fread(fp, [he_w he_h*mul], form, 0, 'ieee-be')';
                imData{cnt-sframe+1}=cast(tmp1,'single');
            end
        else
            for cnt = sframe:lastframe
                %cnt;
                fseek(fp,ofds(cnt),'bof');
                tmp1 = fread(fp, [he_w he_h*mul], form, 0, 'ieee-le')';
                imData{cnt-sframe+1}=cast(tmp1,'single');
            end
        end
        elseif strcmpi(form,'double')
            if(bo)
                for cnt = sframe:lastframe
                    %cnt;
                    fseek(fp,ofds(cnt),'bof');
                    tmp1 = fread(fp, [he_w he_h*mul], form, 0, 'ieee-be.l64')';
                    imData{cnt-sframe+1}=cast(tmp1,'single');
                end
            else
                for cnt = sframe:lastframe
                    %cnt;
                    fseek(fp,ofds(cnt),'bof');
                    tmp1 = fread(fp, [he_w he_h*mul], form, 0, 'ieee-le.l64')';
                    imData{cnt-sframe+1}=cast(tmp1,'single');
                end
            end
    end
            %ieee-le.l64
        
        imData=cell2mat(imData);
        imData=reshape(imData,[he_h*mul,he_w,framenum]);
        fclose(fp);
        display('Finished reading images')
elseif strcmpi(ext,'.hdf5') || strcmpi(ext,'.h5');
    info = hdf5info(path_to_file);
    dims = info.GroupHierarchy.Datasets.Dims;
    if nargin < 2
        sframe = 1;
    end
    if nargin < 3
        num2read = dims(end)-sframe+1;
    end
    num2read = min(num2read,dims(end)-sframe+1);
    imData = h5read(path_to_file,'/mov',[ones(1,length(dims)-1),sframe],[dims(1:end-1),num2read]);
else
    error('Unknown file extension. Only .tiff and .hdf5 files are currently supported');
end