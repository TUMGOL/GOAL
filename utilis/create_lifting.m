%(c) Simon Hawe, Lehrstuhl fuer Datenverarbeitung Technische Universitaet
%Muenchen, 2012. Contact: simon.hawe@tum.de
function K = create_lifting(k,K)
if k==1
    return;
end

if isscalar(K)
    [Vec_pool,~,~] = svd(randn(K));
    K = Vec_pool(1:K,:);
end
n = k*size(K,2);
while 1
    A       = randperm(size(K,1));
    Multi   = randn(1,4);
    Multi   = Multi/norm(Multi);
    New_vec = Multi*K(A(1:4),:);
    New_vec   = New_vec/norm(New_vec);
    if isempty(intersect(New_vec,K,'rows'))
        K = [K;New_vec];
    end
    
    if size(K,1) ==  n
        return;
    end
end
