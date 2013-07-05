function [PSNR,MSE,max_error] = psnr(original, noisy, max_val)
% Calculates PSNR between original image and noisy version
%original = double(uint8(original));
%noisy    = double(uint8(noisy));
if nargin == 2
    max_val = 255;
end

l = 1;
r = l - 1;
original = double(original(l:end-r,l:end-r));
noisy    = double(noisy(l:end-r,l:end-r));
original = double(original);
noisy    = double(noisy);
n=numel(original);

error = abs(original - noisy);

MSE = sum(sum(error.^2))/n;

PSNR = 10*log10(max_val^2/MSE);

max_error = max(abs(error(:)));

fprintf('***********************************\n');
fprintf('Percentage of bad pixels:   %3.4f\n', sum(sum(error>0))/n);
fprintf('Mean Error per Pixel:       %3.4f\n', mean(error(:)));
fprintf('PSNR:                       %3.4f\n', PSNR);
fprintf('MSE:                        %3.4f\n', MSE);
fprintf('RMSE:                       %3.4f\n', sqrt(MSE));
fprintf('Max Error:                  %3.4f\n', max(abs(error(:))));




