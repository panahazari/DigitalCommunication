function [channel] = Channel1(T, sn)
    nT = sn;        
    n = 0:3; 
    h = [1, 1/2, 3/4, -2/7]; 
    figure;
    freqz(h);
    channel = upsample(h,nT);
end