%% Communication Systems Final Project
close all;
clc;
clear all;

qbit = 8;
file_name = '../images/image2.jpg';

[Ztres,r,c,m,n,minval,maxval] = ImagePreProcess_color(file_name, qbit);
title('Original Image');

%%

N = 100000;
[nubmerofGroups, BitStreams] = Convert2Bitstream(N, Ztres, m, n, qbit);


%% Converting binary data to the train of impulses

T = 1;  % Replace with your desired duration between impulses (in seconds)
sn = 32;  % Replace with your desired number of samples between impulses
impulse_train = bit_sequence_to_impulses(BitStreams, sn);

% Plot

% figure;
% stem([0:T/sn:length(impulse_train)/sn - T/sn], impulse_train, 'LineWidth', 2);
% title('Impulse Train with Sample Numbers');
% xlabel('Sample Number');
% ylabel('Amplitude');
% grid on;


%% Modulating impulses with HSPM pulse

T = 1;      % Pulse duration
sn = 32;    % number of samples at each pulse
A1 = 1;     % pulse magnitude
Ts = T/sn; %Every bit has for example 1 sec duration and 32 samples.
pulse = 0:Ts:T - Ts;
hspm_pulse = A1*sin(pi*pulse/T);
hspm_pulse = hspm_pulse./norm(hspm_pulse);
modulated_hspm = conv(hspm_pulse, impulse_train);
modulated_hspm = modulated_hspm(1:length(modulated_hspm)-32);

hspm_pulse_energy = sum(hspm_pulse.^2);
% figure;
% plot(0:T/sn:length(modulated_hspm)/sn - T/sn, modulated_hspm);

%%

%% Modulating impulses with SRRC pulse
alpha = 0.5;
K = 6;
pulse = -K*T:Ts:K*T -Ts;
srrc_pulse = rcosdesign(alpha, 2*K, 32, 'sqrt');
srrc_pulse = srrc_pulse(1:end-1);

srrc_pulse_energy = sum(srrc_pulse.^2);

% figure
% plot(srrc_pulse);
% figure
% freqz(srrc_pulse)

% plot(pulse, srrc_pulse);
modulated_srrc = conv(srrc_pulse, impulse_train);

% figure
% plot(modulated_srrc)

%% Channel:
channel = Channel3();
channel_power_gain = norm(channel);
channel_t = upsample(channel, sn);

figure;
stem(channel);
xlabel('Time (t)');
ylabel('h(t)');
title('Channel impulse response');

figure;
freqz(channel);

modulated_hspm_after_channel = conv(modulated_hspm, channel_t);

modulated_SRRC_after_channel = conv(modulated_srrc, channel_t);

t = 0:1/sn:length(modulated_hspm_after_channel)/sn - 1/sn;


%% Noise
sigma = 0.1;
sigma = [0.005, 0.3, 0];

noise_hspm = sigma(2)*randn(size(modulated_hspm_after_channel));
received_hspm = modulated_hspm_after_channel + noise_hspm;

noise_srrc = sigma(2)*randn(size(modulated_SRRC_after_channel));
received_srrc = modulated_SRRC_after_channel + noise_srrc;

hspm_snr = 10*log10(hspm_pulse_energy./sigma(2));
srrc_snr = 10*log10(srrc_pulse_energy./sigma(2));

%% Matched filter for HSPM
pulse = 0:Ts:T - Ts;
matched_filter = A1*sin(pi*pulse/T);
matched_filter = matched_filter./norm(matched_filter);

matched_output_hspm = conv(received_hspm, matched_filter);

%% Matched filter SRRC
alpha = 0.5;
K = 6;
pulse = -K*T:Ts:K*T - Ts;

srrc_matched_filter = rcosdesign(alpha, 2*K, 32, 'sqrt');
srrc_matched_filter = srrc_matched_filter./norm(srrc_matched_filter);
srrc_matched_filter = srrc_matched_filter(1:end-1);

matched_output_srrc = conv(received_srrc, srrc_matched_filter);

%% Equalizer
h = channel;
eq_zf_hspm = ZF_Equalizer(matched_output_hspm, upsample(h, 32));
eq_zf_srrc = ZF_Equalizer(matched_output_srrc, upsample(h, 32));
eq_mmse_hspm = MMSE_Equalizer(matched_output_hspm, upsample(h, 32), sigma(2));
eq_mmse_srrc = MMSE_Equalizer(matched_output_srrc, upsample(h, 32), sigma(2));


%% Sampling 
sampled_signal_zf_hspm = eq_zf_hspm(sn:sn:length(impulse_train));
sampled_signal_zf_srrc = eq_zf_srrc(2*K*sn:sn:length(impulse_train) + 2*K*sn - sn);

sampled_signal_mmse_hspm = eq_mmse_hspm(sn:sn:length(impulse_train));
sampled_signal_mmse_srrc = eq_mmse_srrc(2*K*sn:sn:length(impulse_train) + 2*K*sn - sn);

%% Detection
detected_zf_hspm = zeros(1, length(sampled_signal_zf_hspm));
detected_zf_hspm(find(sampled_signal_zf_hspm>0))=1;
detected_zf_hspm(find(sampled_signal_zf_hspm<0))=0;


detected_zf_srrc = zeros(1, length(sampled_signal_zf_srrc));
detected_zf_srrc(find(sampled_signal_zf_srrc>0))=1;
detected_zf_srrc(find(sampled_signal_zf_srrc<0))=0;

detected_mmse_hspm = zeros(1, length(sampled_signal_mmse_hspm));
detected_mmse_hspm(find(sampled_signal_mmse_hspm>0))=1;
detected_mmse_hspm(find(sampled_signal_mmse_hspm<0))=0;


detected_mmse_srrc = zeros(1, length(sampled_signal_mmse_srrc));
detected_mmse_srrc(find(sampled_signal_mmse_srrc>0))=1;
detected_mmse_srrc(find(sampled_signal_mmse_srrc<0))=0;



%% Reshaping to matrix
matrix8Bits_zf_hspm = reshape(detected_zf_hspm, 8, []).';
matrix8Bits_zf_srrc = reshape(detected_zf_srrc, 8, []).';
matrix8Bits_mmse_hspm = reshape(detected_mmse_hspm, 8, []).';
matrix8Bits_mmse_srrc = reshape(detected_mmse_srrc, 8, []).';

%% Converting to decimal

decimal_zf_hspm = bi2de(matrix8Bits_zf_hspm);
decimal_zf_srrc = bi2de(matrix8Bits_zf_srrc);
decimal_mmse_hspm = bi2de(matrix8Bits_mmse_hspm);
decimal_mmse_srrc = bi2de(matrix8Bits_mmse_srrc);

%% Reshape the decimal matrix into 8x8 blocks
%% Q15:


SNR_ZF_HSPM = snr(modulated_hspm, eq_zf_hspm(1:1:length(impulse_train)-1))
SNR_MMSE_HSPM = snr(modulated_hspm, eq_mmse_hspm(1:1:length(impulse_train)-1))
SNR_ZF_SRRC = snr(modulated_srrc, eq_zf_srrc(1:1:length(impulse_train) + 2*K*sn-1))
SNR_MMSE_SRRC = snr(modulated_srrc, eq_mmse_srrc(1:1:length(impulse_train) + 2*K*sn-1))

% Recovering ZF HSPM:
resconstrucedZtres = reshape(decimal_zf_hspm, 8, 8, []);
resconstrucedZtres = uint8(resconstrucedZtres);
ImagePostProcess_color(resconstrucedZtres,r,c,m,n,minval,maxval);
title(['Recovered Image using Zero Forcing HSPM, Noise power = ', num2str(sigma(2))]);

% Recovering MMSE HSPM:
resconstrucedZtres = reshape(decimal_mmse_hspm, 8, 8, []);
resconstrucedZtres = uint8(resconstrucedZtres);
ImagePostProcess_color(resconstrucedZtres,r,c,m,n,minval,maxval);
title(['Recovered Image using MMSE HSPM, Noise power = ', num2str(sigma(2))]);

% Recovering ZF SRRC:
resconstrucedZtres = reshape(decimal_zf_srrc, 8, 8, []);
resconstrucedZtres = uint8(resconstrucedZtres);
ImagePostProcess_color(resconstrucedZtres,r,c,m,n,minval,maxval);
title(['Recovered Image using Zero Forcing SRRC, Noise power = ', num2str(sigma(2))]);

% Recovering MMSE SSRC:

resconstrucedZtres = reshape(decimal_mmse_srrc, 8, 8, []);
resconstrucedZtres = uint8(resconstrucedZtres);
ImagePostProcess_color(resconstrucedZtres,r,c,m,n,minval,maxval);
title(['Recovered Image using MMSE SRRC, Noise power = ', num2str(sigma(2))]);










%% 

% %% Q1
% 
% % HSPM
% pulse = 0:Ts:T - Ts;
% hspm_pulse = A1*sin(pi*pulse/T);
% L = length(hspm_pulse);
% figure;
% plot(pulse,hspm_pulse,'LineWidth',1.5)
% xlabel('Time (t)');
% ylabel('g_1(t)');
% title('HSPM Pulse');
% grid on;
% 
% % Plot the magnitude spectrum
% n = 2^nextpow2(L);
% Fs = 1/Ts;
% f = Fs*(-n/2:n/2-1)/n;
% hspm_pulse_freq = fft(hspm_pulse,n);
% hspm_pulse_freq = fftshift(hspm_pulse_freq);
% figure
% freqz(hspm_pulse);
% 
% 
% 
% %% Q1 SRRC
% alpha = 0.5;
% K = 2;
% pulse = -K*T:Ts:K*T - Ts;
% 
% srrc_pulse = rcosdesign(alpha, 2*K, 32, 'sqrt');
% srrc_pulse = srrc_pulse(1:end-1);
% % srrc_pulse = srrc_pulse/norm(srrc_pulse);
% L = length(srrc_pulse);
% figure;
% plot(pulse, srrc_pulse,'LineWidth',1.5)
% xlabel('Time (t)');
% ylabel('g_2(t)');
% xlim([-2 2]);
% title(['SRRC Pulse k = ' int2str(K) , '  a = ', num2str(alpha)]);
% grid on; 
% 
% figure;
% freqz(srrc_pulse);
% 
% 
% %% Q2
% 
% % HSPM Modulated Signal in time for 10 random bits.
% 
% T = 1;      % Pulse duration
% sn = 32;    % number of samples at each pulse
% A1 = 1;     % pulse magnitude
% Ts = T/sn; %Every bit has for example 1 sec duration and 32 samples.
% pulse = 0:Ts:T - Ts;
% hspm_pulse = A1*sin(pi*pulse/T);
% hspm_pulse = hspm_pulse./norm(hspm_pulse);
% modulated_hspm = conv(hspm_pulse, impulse_train);
% modulated_hspm = modulated_hspm(1:length(modulated_hspm)-32);
% figure;
% plot(0:T/sn:length(modulated_hspm)/sn - T/sn, modulated_hspm);
% title('Modulated 10 bits- HSPM');
% xlabel('Time (t)');
% ylabel('s_1(t)');
% grid on;
% 
% %%
% 
% % SRRC Modulated Signal in time for 10 random bits.
% 
% alpha = 0.5;
% K = 6;
% 
% srrc_pulse = rcosdesign(alpha, 2*K, 32, 'sqrt');
% srrc_pulse = srrc_pulse(1:end-1);
% 
% 
% % plot(pulse, srrc_pulse);
% modulated_srrc = conv(srrc_pulse, impulse_train);
% pulse = -K*T:Ts:ceil(length(modulated_srrc)/sn) -K*T - 2*Ts;
% figure
% plot(pulse, modulated_srrc)
% title('Modulated 10 bits- SRRC');
% xlabel('Time (t)');
% ylabel('s_2(t)');
% grid on;
% 
% %% Q3
% % HSPM Spectrum of the modulated signal:
% figure
% freqz(modulated_hspm);
% % SRRC Spectrum of the modulated signal:
% figure
% freqz(modulated_srrc);
% 
% 
% %% Eye Diagrams for HSPM:
% 
% %% Q4 5 7
% 
% % Before convolving with channel:
% fig = eyediagram(modulated_hspm(32*10 + 1:end-32), 32, 1, 16, 'b-');
% set(gca, 'Color','w', 'XColor','black', 'YColor','black');
% set(fig, 'Color', 'white');
% title('Eye diagram for the HSPM modulated signal', 'Color', 'black');
% 
% % After convolving with the channel"
% fig = eyediagram(modulated_hspm_after_channel(32*10 + 1:end-32), 32, 1, 16, 'b-');
% set(gca, 'Color','w', 'XColor','black', 'YColor','black');
% set(fig, 'Color', 'white');
% title('Eye diagram for the HSPM modulated signal with channel', 'Color', 'black');
% 
% %%
% % After adding noise:
% fig = eyediagram(received_hspm(32*10 + 1:end-32), 32, 1, 16, 'b-');
% set(gca, 'Color','w', 'XColor','black', 'YColor','black');
% set(fig, 'Color', 'white');
% title('Eye diagram for the HSPM modulated signal with noise', 'Color', 'black');
% 
% %% Eye Diagrams for SRRC:
% 
% % Before convolving with channel:
% fig = eyediagram(modulated_srrc(32*10 + 1:end-32*10), 32, 32, 0, 'b-');
% set(gca, 'Color','w', 'XColor','black', 'YColor','black');
% set(fig, 'Color', 'white');
% title('Eye diagram for the SRRC modulated signal', 'Color', 'black');
% 
% % After convolving with the channel"
% fig = eyediagram(modulated_SRRC_after_channel(32*10 + 1:end-32*10), 32, 32, 0, 'b-');
% set(gca, 'Color','w', 'XColor','black', 'YColor','black');
% set(fig, 'Color', 'white');
% title('Eye diagram for the SRRC modulated signal with channel', 'Color', 'black');
% %%
% % After adding noise:
% fig = eyediagram(received_srrc(32*10 + 1:end-32*10), 32, 32, 0, 'b-');
% set(gca, 'Color','w', 'XColor','black', 'YColor','black');
% set(fig, 'Color', 'white');
% title('Eye diagram for the SRRC modulated signal with noise', 'Color', 'black');
% 
% %% Matched filter
% 
% %% Q8
% 
% %% HSPM_matched
% 
% pulse = 0:Ts:T - Ts;
% hspm_matched_filter = A1*sin(pi*pulse/T);
% L = length(hspm_matched_filter);
% figure;
% plot(pulse,hspm_matched_filter,'LineWidth',1.5)
% xlabel('Time (t)');
% ylabel('g_1(t)');
% title('Matched filter');
% grid on;
% 
% figure
% freqz(hspm_pulse);
% 
% 
% 
% %% Q1 SRRC
% alpha = 0.5;
% K = 6;
% pulse = -K*T:Ts:K*T - Ts;
% 
% srrc_matched_filter = rcosdesign(alpha, 2*K, 32, 'sqrt');
% srrc_matched_filter = srrc_matched_filter(1:end-1);
% % srrc_pulse = srrc_pulse/norm(srrc_pulse);
% L = length(srrc_matched_filter);
% figure;
% plot(pulse, srrc_matched_filter,'LineWidth',1.5)
% xlabel('Time (t)');
% ylabel('g_2(t)');
% title(['SRRC matched filter k = ' int2str(K) , '  a = ', num2str(alpha)]);
% grid on; 
% 
% figure;
% freqz(srrc_pulse);
% 
% %% Q9
% 
% % 1bit
% fig = eyediagram(matched_output_hspm(32*10 + 1:end-32*10), 32, 32, 0, 'b-');
% set(gca, 'Color','w', 'XColor','black', 'YColor','black');
% set(fig, 'Color', 'white');
% title('Eye diagram for the output of matched filter for HSP modulated signal', 'Color', 'black');
% 
% % 2bit
% fig = eyediagram(matched_output_hspm(32*10 + 1:end-32*10), 64, 32, 0, 'b-');
% set(gca, 'Color','w', 'XColor','black', 'YColor','black');
% set(fig, 'Color', 'white');
% title('Eye diagram for the output of matched filter for HSP modulated signal', 'Color', 'black');
% 
% % 1bit
% fig = eyediagram(matched_output_srrc(32*10 + 1:end-32*10), 32, 32, 0, 'b-');
% set(gca, 'Color','w', 'XColor','black', 'YColor','black');
% set(fig, 'Color', 'white');
% title('Eye diagram for the output of matched filter for SRRC modulated signal', 'Color', 'black');
% 
% % 2bit
% fig = eyediagram(matched_output_srrc(32*10 + 1:end-32*10), 64, 32, 0, 'b-');
% set(gca, 'Color','w', 'XColor','black', 'YColor','black');
% set(fig, 'Color', 'white');
% title('Eye diagram for the output of matched filter for SRRC modulated signal', 'Color', 'black');
% 
% 
% %% Q10
% 
% %% Zero Forcing 
% h = [1, 1/2, 3/4, -2/7];
% impulse = cat(1,1,zeros(800,1));
% impulse_resp = ZF_Equalizer(impulse', h);
% 
% figure;
% stem(impulse_resp);
% 
% figure;
% freqz(impulse_resp);
% 
% %%
% 
% %% Q11 
% 
% h = [1, 1/2, 3/4, -2/7];
% 
% % HSPM
% eq_zf_hspm = ZF_Equalizer(matched_output_hspm, upsample(h, 32));
% fig = eyediagram(eq_zf_hspm(32:floor(length(eq_zf_hspm)/128)), 32, 32, 1, 'b-')
% set(gca, 'Color','w', 'XColor','black', 'YColor','black');
% set(fig, 'Color', 'white');
% title(['Eye diagram for the output of ZF equalizer for HSP modulated signal with noise level = ', num2str(sigma(2))], 'Color', 'black');
% 
% 
% %% SRRC
% 
% eq_zf_srrc = ZF_Equalizer(matched_output_srrc, upsample(h, 32));
% fig = eyediagram(eq_zf_srrc(512:floor(length(eq_zf_srrc)/16)), 32, 32, 0, 'b-');
% set(gca, 'Color','w', 'XColor','black', 'YColor','black');
% set(fig, 'Color', 'white');
% title(['Eye diagram for the output of ZF equalizer for SRRC modulated signal with noise level = ', num2str(sigma(2))], 'Color', 'black');
% 
% %% Q12
% 
% %% MMSE 
% h = [1, 1/2, 3/4, -2/7];
% impulse = cat(1,1,zeros(800,1));
% impulse_resp = MMSE_Equalizer(impulse', h, 0.01);
% 
% figure;
% stem(impulse_resp);
% title('Impulse reponse of the MMSE filter');
% 
% figure;
% freqz(impulse_resp);
% 
% %% Q13
% 
% h = [1, 1/2, 3/4, -2/7];
% 
% % HSPM
% eq_mmse_hspm = MMSE_Equalizer(matched_output_hspm, upsample(h, 32), sigma(2));
% fig = eyediagram(eq_mmse_hspm(32:floor(length(eq_mmse_hspm)/128)), 32, 32, 1, 'b-')
% set(gca, 'Color','w', 'XColor','black', 'YColor','black');
% set(fig, 'Color', 'white');
% title(['Eye diagram for the output of MMSE equalizer for HSP modulated signal with noise level = ', num2str(sigma(2))], 'Color', 'black');
% 
% 
% %% SRRC
% 
% eq_mmse_srrc = MMSE_Equalizer(matched_output_srrc, upsample(h, 32), sigma(2));
% fig = eyediagram(eq_mmse_srrc(512:floor(length(eq_mmse_srrc)/16)), 32, 32, 0, 'b-');
% set(gca, 'Color','w', 'XColor','black', 'YColor','black');
% set(fig, 'Color', 'white');
% title(['Eye diagram for the output of MMSE equalizer filter for SRRC modulated signal with noise level = ', num2str(sigma(2))], 'Color', 'black');
% 
% %%
% 
% %% Q 
% 


