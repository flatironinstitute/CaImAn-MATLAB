function cm = com(A,d1,d2)

% center of mass calculation

nr = size(A,2);
Coor.x = kron(ones(d2,1),(1:d1)');
Coor.y = kron((1:d2)',ones(d1,1));
cm = zeros(nr,2);  % vector for center of mass							   
cm(:,1) = Coor.x'*A(:,1:nr)./sum(A(:,1:nr));
cm(:,2) = Coor.y'*A(:,1:nr)./sum(A(:,1:nr));