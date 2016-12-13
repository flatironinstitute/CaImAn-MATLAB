function AY = mm_fun(A,Y,chunk_size)

% multiply A*Y or A'*Y or Y*A or Y*C' depending on the dimension for loaded
% or memory mapped Y.

memmaped = isobject(Y);
[d1a,d2a] = size(A);

if ~memmaped
    [d1y,d2y] = size(Y);
else
    [d1y,d2y] = size(Y,'Yr');
end

if ~memmaped
    if d1a == d1y
        AY = A'*Y;
    elseif d1a == d2y
        AY = Y*A;
    elseif d2a == d1y
        AY = A*Y;
    elseif d2a == d2y
        AY = Y*A';
    else
        error('matrix dimensions do not agree');
    end
else
    if nargin < 3 || isempty(chunk_size); chunk_size = 2e4; end
    if d1a == d1y
        if nargin < 3 || isempty(chunk_size); chunk_size = 2e4; end
        AY = zeros(d2a,d2y);
        for i = 1:chunk_size:d1a
            AY = AY + A(i:min(i+chunk_size-1,d1a),:)'*double(Y.Yr(i:min(i+chunk_size-1,d1a),:));
        end
    elseif d1a == d2y
        AY = zeros(d1y,d2a);
        for i = 1:chunk_size:d2a
            AY(i:min(i+chunk_size-1,d1),:) = double(Y.Yr(i:min(i+chunk_size-1,d1a),:))*A;
        end
    elseif d2a == d1y
        if nargin < 3 || isempty(chunk_size); chunk_size = 2e4; end
        AY = zeros(d1a,d2y);
        for i = 1:chunk_size:d2a
            AY = AY + A(:,i:min(i+chunk_size-1,d1a))'*double(Y.Yr(i:min(i+chunk_size-1,d1a),:));
        end
    elseif d2a == d2y
        AY = zeros(d1y,d1a);
        At = A';
        for i = 1:chunk_size:d1a
            AY(i:min(i+chunk_size-1,d1),:) = double(Y.Yr(i:min(i+chunk_size-1,d1a),:))*At;
        end
    else
        error('matrix dimensions do not agree');
    end
end