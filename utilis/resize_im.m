function x = resize_im(x,newsiz)

siz = size(x);
N = prod(siz);
M = prod(newsiz);

% DCT transform
x = dctn(x);

% Crop the DCT coefficients
if M < N
    x = x(1:newsiz(1),1:newsiz(2));
    mul = sqrt(M/N);
else
    x = padarray(x,newsiz-siz,0,'post');
    mul = sqrt(N/M);
end
% inverse DCT transform
x = idctn(x)*mul;


end
