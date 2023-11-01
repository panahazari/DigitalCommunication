function [modulatedBitStream, t] = HSPM(bitStream, T, sn, A1) %Half Sine Pulse Modulation
    Ts = T/sn; %Every bit has for example 1 sec duration and 32 samples.
    t = Ts:Ts:T*length(bitStream);    
    pulse = 0:Ts:1 - Ts;
     
    for k = 1 : length(bitStream)
        if bitStream(k) == 1
            modulatedBitStream((k-1)*sn +1 : k*sn) =  A1*sin(pi*pulse/T);
        else
            modulatedBitStream((k-1)*sn +1 : k*sn) =  -A1*sin(pi*pulse/T);
        end
    end

end