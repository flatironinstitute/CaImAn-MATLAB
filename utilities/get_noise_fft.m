 function [sn,psdx,ff] = get_noise_fft(Y,options)
        
        defoptions = CNMFSetParms;
        if nargin < 2 || isempty(options); options = defoptions; end

        if ~isfield(options,'noise_range'); options.noise_range = defoptions.noise_range; end
        range_ff = options.noise_range;
        if ~isfield(options,'noise_method'); options.noise_method = defoptions.noise_method; end
        method = options.noise_method;
        if ~isfield(options,'block_size'); options.block_size = defoptions.block_size; end                
        block_size = options.block_size;
        if ~isfield(options,'split_data'); options.split_data = defoptions.split_data; end
        split_data = options.split_data;
        if ~isfield(options,'max_timesteps') || isempty(options.max_timesteps); 
            options.max_timesteps = defoptions.max_timesteps;
        end
        
        dims = ndims(Y);
        sizY = size(Y);
        N = min(sizY(end),options.max_timesteps);
        if N < sizY(end)
           %Y = reshape(Y,prod(sizY(1:end-1)),[]);
           switch ndims(Y), 
               case 3, 
                    Y(:,:,N+1:end) = [];
               case 4, 
                    Y(:,:,:,N+1:end) = [];
           end
        end
        
        Fs = 1;        
        ff = 0:Fs/N:Fs/2;
        indf=ff>range_ff(1);
        indf(ff>range_ff(2))=0;
        if dims > 1
            d = prod(sizY(1:dims-1));
            Y = reshape(Y,d,N);
            Nb = prod(block_size);
            SN = cell(ceil(d/Nb),1);
            PSDX = cell(ceil(d/Nb),1);
            if ~split_data
                for ind = 1:ceil(d/Nb); 
                    xdft = fft(Y((ind-1)*Nb+1:min(ind*Nb,d),:),[],2); 
                    xdft = xdft(:,1: floor(N/2)+1); % FN: floor added.
                    psdx = (1/(Fs*N)) * abs(xdft).^2;
                    psdx(:,2:end-1) = 2*psdx(:,2:end-1) + eps;
                    %SN{ind} = mean_psd(psdx(:,indf),method);
                    switch lower(method)
                        case 'mean'
                            SN{ind}=sqrt(mean(psdx(:,indf)/2,2));
                        case 'median'
                            SN{ind}=sqrt(median(psdx(:,indf)/2,2));
                        case 'logmexp'
                            SN{ind} = sqrt(exp(mean(log(psdx(:,indf)/2),2)));
                        otherwise
                            error('unknown method for averaging noise..')
                    end
                    PSDX{ind} = psdx;
                end
            else
                nc = ceil(d/Nb);
                Yc = mat2cell(Y,[Nb*ones(nc-1,1);d-(nc-1)*Nb],N);
                parfor ind = 1:ceil(d/Nb); 
                    xdft = fft(Yc{ind},[],2); 
                    xdft = xdft(:,1:floor(N/2)+1);
                    psdx = (1/(Fs*N)) * abs(xdft).^2;
                    psdx(:,2:end-1) = 2*psdx(:,2:end-1) + eps;
                    Yc{ind} = [];
                    switch lower(method)
                        case 'mean'
                            SN{ind}=sqrt(mean(psdx(:,indf)/2,2));
                        case 'median'
                            SN{ind}=sqrt(median(psdx(:,indf)/2,2));
                        case 'logmexp'
                            SN{ind} = sqrt(exp(mean(log(psdx(:,indf)/2),2)));
                        otherwise
                            error('unknown method for averaging noise..')
                    end
                    
                end
            end
            sn = cell2mat(SN);
        else
            xdft = fft(Y);
            xdft = xdft(:,1:floor(N/2)+1);
            psdx = (1/(Fs*N)) * abs(xdft).^2;
            psdx(:,2:end-1) = 2*psdx(:,2:end-1) + eps;
            switch lower(method)
                case 'mean'
                    sn = sqrt(mean(psdx(:,indf)/2,2));
                case 'median'
                    sn = sqrt(median(psdx(:,indf)/2,2));
                case 'logmexp'
                    sn = sqrt(exp(mean(log(psdx(:,indf)/2),2)));
                otherwise
                    error('unknown method for averaging noise..')
            end
        end
        psdx = cell2mat(PSDX);
        if dims > 2
            sn = reshape(sn,sizY(1:dims-1));
        end
 end
