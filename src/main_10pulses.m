%% Communication Systems Final Project
close all;
clc;
clear all;

qbit = 8;

%% Converting binary data to the train of impulses

T = 1;  % Replace with your desired duration between impulses (in seconds)
sn = 32;  % Replace with your desired number of samples between impulses

%BitStreams = [ 1 1 0 1 0 1 1 0 0];
BitStreams = [1];
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
channel_t = Channel1(T, 32);
t = 0:1/sn:length(channel_t)/sn - 1/sn;

% figure;
% stem(t, channel_t);
% xlabel('Time (t)');
% ylabel('h(t)');
% title('Channel impulse response');

% figure;
% freqz(channel_t);

modulated_hspm_after_channel = conv(modulated_hspm, channel_t);

modulated_SRRC_after_channel = conv(modulated_srrc, channel_t);

t = 0:1/sn:length(modulated_hspm_after_channel)/sn - 1/sn;


%% Noise
sigma = 0.1;
sigma = [0.005, 0, 0];

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
h = [1, 1/2, 3/4, -2/7];
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

figure;
t = -1:1/sn:(1/sn)*length(eq_mmse_hspm)-1 -1/sn;
plot(t, eq_mmse_hspm);
grid on;
xlim([-2 2]);
title("The output of equalizer-HSP");

figure; 
t = -2*K:1/sn:length(eq_mmse_srrc)*(1/sn)-2*K - 1/sn;
plot(t,eq_mmse_srrc);
xlim([-6 6]);
grid on;
title("The output of equalizer-SRRC");
