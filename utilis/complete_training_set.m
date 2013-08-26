%(c) Simon Hawe, Lehrstuhl fuer Datenverarbeitung Technische Universitaet
%Muenchen, 2012. Contact: simon.hawe@tum.de
%Images    => Cell of images used to extract patches or path to folder
%dim       => dimension of single patch
%n_patches => Total Number of patches
function [S] = complete_training_set(Images, n_patches, Patch_size)
S = [];
if  ischar(Images)
    Image_list = dir(Images);
    Image_names = {};
    k = 0;
    for i=1:numel(Image_list)
        if ~Image_list(i).isdir && ~strcmp(Image_list(i).name,'Thumbs.db')
            k = k + 1;
            Image_names{k} = [Images,Image_list(i).name];
        end
    end
    n_patches = ceil(n_patches/k);
    for i = 1:k
        Im = imread(Image_names{i});
        if ndims(Im) == 3
            Im = rgb2gray(Im);
        end
        Im = double(Im);
        Sc = extract_training_set(Im, n_patches, Patch_size);
        S  = [S,Sc];
        %S  = [S,extract_training_set(Im, n_patches, Patch_size, method)];
    end
else
    n_patches = ceil(n_patches/numel(Images));
    for i = 1:numel(Images)
        S = [S,extract_training_set(Images{i}, n_patches, Patch_size)];
    end
end



