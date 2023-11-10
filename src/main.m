%% Communication Systems Final Project
close all;
clc;
clear all;

qbit = 8;
file_name = '../images/image1.png';

[Ztres,r,c,m,n,minval,maxval] = ImagePreProcess_gray(file_name, qbit);


%%

N = 10;
[nubmerofGroups, BitStreams] = Convert2Bitstream(N, Ztres, m, n, qbit);


%% Converting binary data to the train of impulses

T = 1;  % Replace with your desired duration between impulses (in seconds)
sn = 32;  % Replace with your desired number of samples between impulses
impulse_train = bit_sequence_to_impulses(BitStreams, sn);
%impulse_train = upsample([1 0 0 1 1], 32);
figure;
stem([0:T/sn:length(impulse_train)/sn - T/sn], impulse_train, 'LineWidth', 2);
title('Impulse Train with Sample Numbers');
xlabel('Sample Number');
ylabel('Amplitude');
grid on;


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
figure;
plot(0:T/sn:length(modulated_hspm)/sn - T/sn, modulated_hspm);

%%

%% Modulating impulses with SRRC pulse
alpha = 0.5;
K = 4;
pulse = -K*T:Ts:K*T -Ts;
srrc_pulse = rcosdesign(alpha, 2*K, 32, 'sqrt');
srrc_pulse = srrc_pulse(1:end-1);

figure
plot(srrc_pulse);
figure
freqz(srrc_pulse)

% plot(pulse, srrc_pulse);
modulated_srrc = conv(srrc_pulse, impulse_train);

figure
plot(modulated_srrc)

%% Channel:
channel_t = Channel1(T, 32);
t = 0:1/sn:length(channel_t)/sn - 1/sn;
figure;
stem(t, channel_t);

%freqz(channel_t);

modulated_hspm_after_channel = conv(modulated_hspm, channel_t);

modulated_SRRC_after_channel = conv(modulated_srrc, channel_t);

t = 0:1/sn:length(modulated_hspm_after_channel)/sn - 1/sn;
%plot(t, modulated_hspm_after_channel);

%% Noise
sigma = 0.1;
sigma = [0.005, 0.1, 0.5];

noise_hspm = sigma(1)*randn(size(modulated_hspm_after_channel));
received_hspm = modulated_hspm_after_channel + noise_hspm;

noise_srrc = sigma(1)*randn(size(modulated_SRRC_after_channel));
received_srrc = modulated_SRRC_after_channel + noise_srrc;

%% Eye Diagrams for HSPM:

% Before convolving with channel:
eyediagram(modulated_hspm, 32, 1, 16);
title('Eye diagram for the modulated signal');

% After convolving with the channel"
eyediagram(modulated_hspm_after_channel, 32, 1, 16);
title('Eye diagram for the modulated signal with channel');

% After adding noise:
eyediagram(received_hspm, 32, 1, 16);
title('Eye diagram for the modulated signal with noise');

%% Eye Diagrams for SRRC:

% Before convolving with channel:
eyediagram(modulated_srrc(32*10:end-32*10), 2*32);
title('Eye diagram for the modulated signal');

% After convolving with the channel"
eyediagram(modulated_SRRC_after_channel(32*10:end-32*10), 2*32);
title('Eye diagram for the modulated signal with channel');

% After adding noise:
eyediagram(received_srrc(32*10:end-32*10), 2*32);
title('Eye diagram for the modulated signal with noise');

%% Q1

% HSPM
pulse = 0:Ts:T - Ts;
hspm_pulse = A1*sin(pi*pulse/T);
L = length(hspm_pulse);
figure;
plot(pulse,hspm_pulse,'LineWidth',1.5)
xlabel('Time (t)');
ylabel('g_1(t)');
title('HSPM Pulse');
grid on;

% Plot the magnitude spectrum
n = 2^nextpow2(L);
Fs = 1/Ts;
f = Fs*(-n/2:n/2-1)/n;
hspm_pulse_freq = fft(hspm_pulse,n);
hspm_pulse_freq = fftshift(hspm_pulse_freq);
figure
freqz(hspm_pulse);



%% Q1 SRRC
alpha = 0.5;
K = 6;
pulse = -K*T:Ts:K*T - Ts;

srrc_pulse = 10*rcosdesign(alpha, 2*K, 32, 'sqrt');
srrc_pulse = srrc_pulse(1:end-1);
% srrc_pulse = srrc_pulse/norm(srrc_pulse);
L = length(srrc_pulse);
figure;
plot(pulse, srrc_pulse,'LineWidth',1.5)
xlabel('Time (t)')
ylabel('g_2(t)')
title(['SRRC Pulse k = ' int2str(K) , '  a = ', num2str(alpha)]);
grid on; 

figure;
freqz(srrc_pulse);


%% Q2

% HSPM Modulated Signal in time for 10 random bits.






% SRRC Modulated Signal in time for 10 random bits.


%% Q3

% HSPM Spectrum of the modulated signal:

% SRRC Spectrum of the modulated signal:

%% Q4:

% Eye Diagram for HSPM:

% Eye Diagram for SRRC:

%% Q5

% Frequency response of the channel

% impulse response of the channel




%% Q6:

% HSPM Channel Output Eye Diagram


 

