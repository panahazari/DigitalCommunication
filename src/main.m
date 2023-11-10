%% Communication Systems Final Project
close all;
clc;
clear all;

qbit = 8;
file_name = '../images/image2.jpg';

[Ztres,r,c,m,n,minval,maxval] = ImagePreProcess_gray(file_name, qbit);

N = 10;
[nubmerofGroups, BitStreams] = Convert2Bitstream(N, Ztres, m, n, qbit);

T = 1;      % Pulse duration
sn = 32;    % number of samples at each pulse
A1 = 1;     % pulse magnitude
testData = [ 1 0 1 0 0 1 1];
[modulatedSignal, t] = HSPM(testData, T, sn, A1);

%plot(t, modulatedSignal);

alpha = 0.9;
K = 6;
Ts = T/sn; %Every bit has for example 1 sec duration and 32 samples.
%t = Ts:Ts:T*length(bitStream);    
% pulse = -K*T:Ts:K*T - Ts;
% 
% 
% 
s = rcosdesign(alpha, 2*K, 32);

% plot(s)

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

% Plot the magnitude spectrum
n = 2^nextpow2(L);
Fs = 1/Ts;
f = Fs*(-n/2:n/2-1)/n;
hspm_pulse_freq = fft(hspm_pulse,n);
hspm_pulse_freq = fftshift(hspm_pulse_freq);
figure
plot(f,10*log10(abs(hspm_pulse_freq)),'LineWidth',1.5);
xlabel({'Frequency (Hz)'},'Interpreter','latex','FontSize',14);
ylabel({'$|G_1(f)|$'},'Interpreter','latex','FontSize',14)
title('Magnitude of the Frequency reponse of HSPM pulse', 'Interpreter', 'latex','FontSize',14);
xlim([-30 30])

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
xlim([-30 30])

% SRRC
pulse = 0:Ts:T - Ts;
srrc_pulse = rcosdesign(alpha, 2*K, 32);
srrc_pulse = srrc_pulse./norm(srrc_pulse);
L = length(srrc_pulse);
figure;
plot(srrc_pulse,'LineWidth',1.5)
xlabel({'Time $(t)$'},'Interpreter','latex','FontSize',14)
ylabel({'$g_2(t)$'},'Interpreter','latex','FontSize',14)
title(['SRRC Pulse k = ' int2str(K) , '$\alpha$ = ', num2str(alpha)],'Interpreter','latex','FontSize',14);

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
title(['Magnitude Frequency reponse of SRRC pulse  k = ' int2str(K) , '$\alpha$ = ', num2str(alpha)], 'Interpreter', 'latex','FontSize',14);
xlim([-20 20])

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
title(['Angle Frequency reponse of SRRC pulse  k = ' int2str(K) , '$\alpha$ = ', num2str(alpha)], 'Interpreter', 'latex','FontSize',14);
xlim([-20 20])


