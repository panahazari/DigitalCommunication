function [output] = MMSE_Equalizer(signal, channel, sigma)

N = 2^nextpow2(length(signal));

E_b = 1; %pulse energy: since we have normalized the pulse enery is equal to 1



channel_f = fft(channel, N);
signal_f = fft(signal, N);

denom = (abs(channel_f)).^2 + sigma/E_b;
mmse_eq = conj(channel_f)./(denom);


output_f = mmse_eq.*(signal_f);

output = ifft(output_f, N);


