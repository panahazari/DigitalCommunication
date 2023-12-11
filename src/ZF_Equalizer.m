function [output] = ZF_Equalizer(signal, channel)

    N = 2^nextpow2(length(signal));
    
    channel_f = fft(channel, N);
    signal_f = fft(signal, N);
    
    zf_eq = 1./channel_f;
    
    output_f = zf_eq.*(signal_f);
    
    output = ifft(output_f, N);

end
