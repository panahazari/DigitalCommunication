%% Communication Systems Final Project
close all;
clc;
clear all;

qbit = 8;
file_name = '../images/image2.jpg';

[Ztres,r,c,m,n,minval,maxval] = ImagePreProcess_gray(file_name, qbit);


%%

N = 1;
[nubmerofGroups, BitStreams] = Convert2Bitstream(N, Ztres, m, n, qbit);


%% Converting binary data to the train of impulses

T = 1;  % Replace with your desired duration between impulses (in seconds)
sn = 32;  % Replace with your desired number of samples between impulses
impulse_train = bit_sequence_to_impulses(BitStreams(1, 1:100), sn);
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
% alpha = 0.5;
% K = 6;
% pulse = -K*T:Ts:K*T;
% srrc_pulse = rcosdesign(alpha, 2*K, 32, 'sqrt');
% srrc_pulse = srrc_pulse./norm(srrc_pulse);
% srrc(pulse) = srrc_pulse;
% srrc_pulse = sigshift(srrc_pulse, k);
% % plot(pulse, srrc_pulse);
% lated_srrc = conv(srrc_pulse, impulse_train);
% figure;
% plot(-K:T/sn:length(modulated_srrc)/sn - T/sn, modulated_srrc);


%% Channel:
channel_t = Channel1(T, 32);
t = 0:1/sn:length(channel_t)/sn - 1/sn;
figure;
stem(t, channel_t);

%freqz(channel_t);

modulated_hspm_after_channel = conv(modulated_hspm, channel_t);

%modulated_SRRC_after_channel = conv(modulated_srrc)

t = 0:1/sn:length(modulated_hspm_after_channel)/sn - 1/sn;
%plot(t, modulated_hspm_after_channel);

%% Noise
sigma = 0.1;
sigma = [0.01, 0.1, 0.5];

noise_hspm = sigma(1)*randn(size(modulated_hspm_after_channel));
received_hspm = modulated_hspm_after_channel + noise_hspm;

% noise_srrc = sigma*randn(size(modulated_SRRC_after_channel));
% received_srrc = modulated_SRRC_after_channel + noise_hspm;

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

%% Q1

% HSPM
pulse = 0:Ts:T - Ts;
hspm_pulse = A1*sin(pi*pulse/T);
L = length(hspm_pulse);
figure;
plot(pulse,hspm_pulse,'LineWidth',1.5)
xlabel({'Time $(t)$'},'Interpreter','latex','FontSize',14)
ylabel({'$g_1(t)$'},'Interpreter','latex','FontSize',14)
title('HSPM Pulse','Interpreter','latex','FontSize',14);
grid on;

% Plot the magnitude spectrum
n = 2^nextpow2(L);
Fs = 1/Ts;
f = Fs*(-n/2:n/2-1)/n;
hspm_pulse_freq = fft(hspm_pulse,n);
hspm_pulse_freq = fftshift(hspm_pulse_freq);
figure
freqz(hspm_pulse)
title('Freq. resp. of SSRC Pulse')

plot(f,10*log10(abs(hspm_pulse_freq)),'LineWidth',1.5);
xlabel({'Frequency (Hz)'},'Interpreter','latex','FontSize',14);
ylabel({'$|G_1(f)|$'},'Interpreter','latex','FontSize',14)
title('Magnitude of the Frequency reponse of HSPM pulse', 'Interpreter', 'latex','FontSize',14);
%xlim([-30 30])
grid on;

% Plot the angle spectrum
n = 2^nextpow2(L);
Fs = 1/Ts;
f = Fs*(-n/2:n/2-1)/n;
hspm_pulse_freq = fft(hspm_pulse,n);
hspm_pulse_freq = fftshift(hspm_pulse_freq);
figure
plot(f,angle(hspm_pulse_freq),'LineWidth',1.5);
xlabel({'Frequency (Hz)'},'Interpreter','latex','FontSize',14);
ylabel({'$\angle G_1(f)$'},'Interpreter','latex','FontSize',14)
title('Angle of the Frequency reponse of HSPM pulse', 'Interpreter', 'latex','FontSize',14);
%xlim([-30 30])
grid on;


%% SRRC
alpha = 0.5;
K = 6;
pulse = -K*T:Ts:K*T;

srrc_pulse = rcosdesign(alpha, 2*K, 32, 'sqrt');
srrc_pulse = srrc_pulse/norm(srrc_pulse);
L = length(srrc_pulse);
figure;
plot(pulse, srrc_pulse,'LineWidth',1.5)
xlabel({'Time $(t)$'},'Interpreter','latex','FontSize',14)
ylabel({'$g_2(t)$'},'Interpreter','latex','FontSize',14)
title(['SRRC Pulse k = ' int2str(K) , '  $\alpha$ = ', num2str(alpha)],'Interpreter','latex','FontSize',14);
grid on;

% Plot the magnitude spectrum
n = 2^nextpow2(L);
Fs = 1/Ts;
f = Fs*(-n/2:n/2-1)/n;
srrc_pulse_freq = fft(srrc_pulse,n);
srrc_pulse_freq = fftshift(srrc_pulse_freq);
figure
plot(f,10*log10(abs(srrc_pulse_freq)),'LineWidth',1.5);
xlabel({'Frequency (Hz)'},'Interpreter','latex','FontSize',14);
ylabel({'$|G_2(f)|$'},'Interpreter','latex','FontSize',14)
title(['Magnitude Frequency reponse of SRRC pulse  k = ' int2str(K) , '  $\alpha$ = ', num2str(alpha)], 'Interpreter', 'latex','FontSize',14);
%xlim([-20 20])
grid on;

% Plot the angle spectrum
n = 2^nextpow2(L);
Fs = 1/Ts;
f = Fs*(-n/2:n/2-1)/n;
srrc_pulse_freq = fft(srrc_pulse,n);
srrc_pulse_freq = fftshift(srrc_pulse_freq);
figure
plot(f,angle(srrc_pulse_freq),'LineWidth',1.5);
xlabel({'Frequency (Hz)'},'Interpreter','latex','FontSize',14);
ylabel({'$\angle G_2(f)$'},'Interpreter','latex','FontSize',14)
title(['Angle Frequency reponse of SRRC pulse  k = ' int2str(K) , '  $\alpha$ = ', num2str(alpha)], 'Interpreter', 'latex','FontSize',14);
%xlim([-20 20])
grid on;

%% Q2

% HSPM Modulated Signal in time for 10 random bits.






% SRRC Modulated Signal in time for 10 random bits.


%% Q3

% HSPM Spectrum of the modulated signal:

% SRRC Spectrum of the modulated signal:

%% Q4:

% Eye Diagram for HSPM:

% Eye Diagram for SRRC:


 

