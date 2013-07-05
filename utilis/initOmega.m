function Omega = initOmega(X,rows)
Omega = zeros(rows,size(X,1));
d = size(X,1)-2;
for i=1:rows
    sel = randperm(length(X));
    Omega(i,:) = null([X(:,sel(1:d))';ones(1,size(X,1))])';
end
% [U,S,V]=svd(Omega);
% S(S~=0) = 1;
% Omega =U*S*V';