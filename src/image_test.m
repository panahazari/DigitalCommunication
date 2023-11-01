close all;
clear all;


qbits = 8	 %% valid options are 8 and 16

filename = 'onion.png';


% Show a colored image in color

[Ztres,r,c,m,n,minval,maxval] = ImagePreProcess_color(filename,qbits);

ImagePostProcess_color(Ztres,r,c,m,n,minval,maxval);


% Show a colored image as 3 separate (gray) layers

[Ztres,r,c,m,n,minval,maxval] = ImagePreProcess_gray(filename,qbits);

ImagePostProcess_gray(Ztres,r,c,m,n,minval,maxval);

