Folder = '..\..\mxTV\Results\';
%Folder = 'C:\BM3D\';
Operators = dir(Folder);
%fileID = fopen('BM3D_PSNR_SSIM.txt','w');
%fileID = fopen('FoE_PSNR_SSIM.txt','w');
%fileID = fopen('KSVD_PSNR_SSIM.txt','w');
fileID = fopen('TV_PSNR_SSIM.txt','w');
for i=1:numel(Operators)
    if ~Operators(i).isdir && ~strcmp(Operators(i).name,'Thumbs.db')
        Crnt_img = (imread([Folder,Operators(i).name]));
        pos = strfind(Operators(i).name,'_');
        if isempty(pos)
			 fprintf(fileID,'**********************************\n');
			 fprintf(fileID,'**********************************\n');
             GT = Crnt_img;
			 continue;
        end
        pos = strfind(Operators(i).name,'.');
		 
		SIDX=ssim_index(double(uint8(Crnt_img)),double(uint8(GT)));
        PSNR = psnr(double(uint8(Crnt_img)),double(uint8(GT)),255);
		 
		fprintf(fileID,'SSIM %s \t \t%.3f ',Operators(i).name(1:pos(1)-1),SIDX);
		fprintf(fileID,'PSNR \t \t %.2f\n',PSNR);
    end
end
fclose(fileID);