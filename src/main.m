%% Communication Systems Final Project
close all;
clc;
clear all;

qbit = 8;
file_name = '../images/image1.png';

[Ztres,r,c,m,n,minval,maxval] = ImagePreProcess_gray(file_name, qbit);

N = 1;
[nubmerofGroups,BitStreams] = Convert2Bitstream(N, Ztres, m, n, qbit);

T = 1;      % Pulse duration
sn = 32;    % number of samples at each pulse
A1 = 1;     % pulse magnitude
[modulatedSignal, t] = HSPM(BitStreams(1,:), T, sn, A1);

