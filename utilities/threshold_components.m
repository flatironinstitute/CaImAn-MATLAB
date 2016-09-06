function Ath = threshold_components(A,options)

% post processing of spatial components
% for each component perform the following:
%   (i)     perform median filtering 
%   (ii)    threshold components
%   (iii)   perform morphological closing
%   (iv)    extract largest connected component

% The thresholding can be determined in two ways:
% options.thr_method = 'max': For every component pixel values below
%       options.thr_method are discarded
% options.thr_method = 'nrg': The boundary is drawn around the value above
%        which a specified fraction of energy is explained (default 99%)

% Written by:
% Eftychios A. Pnevmatikakis, Simons Foundation, 2015

    defoptions = CNMFSetParms;
    if nargin < 2 || isempty(options); options = defoptions; end
%     defoptions.nrgthr = 0.9999;              % energy threshold
%     defoptions.clos_op = strel('square',3);  % morphological operator for closing
%     defoptions.medw = [3,3];                 % size of median filter
    
    if ~isfield(options,'thr_method') || isempty(options.thr_method); options.thr_method = defoptions.thr_method; end
    if ~isfield(options,'nrgthr') || isempty(options.nrgthr); options.nrgthr = defoptions.nrgthr; end
    if ~isfield(options,'maxthr') || isempty(options.maxthr); options.maxthr = defoptions.maxthr; end
    if ~isfield(options,'clos_op') || isempty(options.clos_op); options.clos_op = defoptions.clos_op; end
    if ~isfield(options,'conn_comp') || isempty(options.conn_comp); options.conn_comp = defoptions.conn_comp; end % extract largest connected component
    if ~isfield(options,'medw') || isempty(options.medw); options.medw = defoptions.medw; end
    if ~isfield(options,'d3') || isempty(options.d3); options.d3 = 1; end
    
    if options.d3 > 1; options.thr_method = 'nrg'; end
    [d,nr] = size(A);
    Ath = spalloc(d,nr,nnz(A));
    Ath(:,nr-options.nb+1:nr) = A(:,nr-options.nb+1:nr);
    indf = cell(nr,1);
    valf = cell(nr,1);
    parfor i = 1:nr-options.nb
        A_temp = reshape(full(A(:,i)),options.d1,options.d2,options.d3);
        for z = 1:options.d3
            A_temp(:,:,z) = medfilt2(A_temp(:,:,z),options.medw);
        end
        if strcmpi(options.thr_method,'nrg')
            A_temp = A_temp(:);
            [temp,ind] = sort(A_temp(:).^2,'ascend'); 
            temp =  cumsum(temp);
            ff = find(temp > (1-options.nrgthr)*temp(end),1,'first');
            BW = zeros(options.d1,options.d2,options.d3);
            BW(ind(ff:d)) = 1;
            for z = 1:options.d3
                BW(:,:,z) = imclose(BW(:,:,z),options.clos_op);
            end
            [L,NUM] = bwlabeln(BW,8*(options.d3==1) + 6*(options.d3~=1));
            if NUM > 0
                nrg = zeros(NUM,1);
                for l = 1:NUM
                    ff = (L==l);
                    nrg(l) = sum(A_temp(ff).^2);
                end
                [~,indm] = max(nrg);
                ff = find(L==indm);
                indf{i} = ff;
                valf{i} = A_temp(ff);
            else
                valf{i} = 0;
            end
        else
            if ~strcmpi(options.thr_method,'max')
                fprintf('Unknown thresholding method, using thresholding based on max value. \n');
            end
            A_temp(A_temp<options.maxthr*max(A_temp(:))) = 0;
            BW = imclose(A_temp>0,options.clos_op);
            if options.conn_comp;
                BW = bwareafilt(BW,1);            
            end
            ff = find(BW>0);            
            indf{i} = ff;
            if ~isempty(ff)
                valf{i} = A_temp(ff);
            else
                valf{i} = 0;
            end            
        end
    end   
    for i = 1:nr-options.nb
        if ~isempty(indf{i})
            Ath(indf{i},i) = valf{i};
        end
    end
end