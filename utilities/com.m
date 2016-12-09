function cm = com(A,d1,d2,d3)

% center of mass calculation
% inputs:
% A: d X nr matrix, each column in the spatial footprint of a neuron
% d1, d2, d3: the dimensions of the 2-d (or 3-d) field of view

% output:
% cm: nr x 2 (or 3) matrix, with the center of mass coordinates

if nargin < 4
    d3 = 1;
end
if d3 == 1
    ndim = 2;
else
    ndim = 3;
end

nr = size(A,2);
Coor.x = kron(ones(d2*d3,1),double(1:d1)');
Coor.y = kron(ones(d3,1),kron(double(1:d2)',ones(d1,1)));
Coor.z = kron(double(1:d3)',ones(d2*d1,1));
cm = [Coor.x, Coor.y, Coor.z]'*A/spdiags(sum(A)',0,nr,nr);
cm = cm(1:ndim,:)';