function [channel] = Channel2(T, sn)
    nT = sn;        
    n = 0:3; 
    h = zeros(1,26);
    h(1) = 0.5;
    h(2) = 1;
    h(4) = 0.63;
    h(9) = 0.25;
    h(13) = 0.16;
    h(26) = 0.1;
    figure;
    freqz(h);
    channel = upsample(h,nT);
end
