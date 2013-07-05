function TV_Mat = create_tvmat(sz,~)

Tv_mat                = eye(sz);
Tv_mat(sz+1:sz+1:end) = -1;
%Tv_mat(end,:) = [];

Tv_y = kron(eye(sz),Tv_mat);
Tv_x = eye(sz^2) - diag(ones(sz^2-sz,1),sz);
%Tv_x(end-sz+1:end,:)=[];

TV_Mat = [Tv_x;Tv_y];
return
z = sum(TV_Mat ~= 0,2);
TV_Mat(z==1,:)=[];
TV_Mat = bsxfun(@times,TV_Mat,1./sqrt(sum(TV_Mat.^2,2)));
if nargin == 2
    a = kron(ones(1,sz-1),[ones(1,sz-1),0]);
    a = a(1:end-1);
    Tv_xy=eye(sz^2) - diag(a,sz+1);
    TV_Mat = [TV_Mat;Tv_xy];
end
