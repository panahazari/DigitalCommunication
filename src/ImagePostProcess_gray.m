function []=ImagePostProcess_gray(Ztres,r,c,m,n,minval,maxval)

%% invert the reshaping operation
newZt = reshape(permute(reshape(Ztres,8,8,r,c), [1 3 2 4]), m,n);

%%%%%%%%%%%%%% IMAGE POST-PROCESSING %%%%%%%%%%%%%%%%
temp=im2double(newZt)*(maxval-minval)+minval;
fun=@idct2;
newZ=blkproc(temp,[8 8],fun);

figure('Renderer', 'painters', 'Position', [10 10 1000 1000]);
imshow(newZ);
