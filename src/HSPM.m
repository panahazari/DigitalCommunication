function [modulatedBitStream, t] = HSPM(bitStream, T, sn, A1) %Half Sine Pulse Modulation
    Ts = 1/sn; %Every bit has for example 1 sec duration and 32 samples.
    t = Ts:Ts:T*length(bitStream);    
    pulse = 0:Ts:T - Ts;
     
    for k = 1 : length(bitStream)
        if bitStream(k) == 1
            modulatedBitStream((k-1)*sn*T +1 : k*sn*T) =  A1*sin(pi*pulse/T);
        else
            modulatedBitStream((k-1)*sn*T +1 : k*sn*T) =  -A1*sin(pi*pulse/T);
        end

    end
    
    modulatedBitStream = modulatedBitStream./norm(modulatedBitStream);
end