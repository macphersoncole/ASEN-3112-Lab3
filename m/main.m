%% ASEN 3112 Structures Lab #3 - main.m
%
%   Author: Andrew Thompson
%   Created: 11/22/20 Edited: 11/22/20

clc; clear all; close all;

%% Import & Clean Data
rawData = importdata("Lab3Data.txt");
data = rawData.data;
time = data(:, 1) - data(1, 1);
accCh0  = data(:,2);
accCh1  = data(:,3);
accCh2  = data(:,4);
accCh3  = data(:,5);
dispCh0 = data(:,6);
dispCh1 = data(:,7);
dispCh2 = data(:,8);
dispCh3 = data(:,9);
vibro   = data(:,10);

%% Excitation Frequency
% Sampling frequency 
Fs = 2500; 
% Sampling Period
T = 1/Fs;
% Sample length
L = length(time);

idx3 = 1;
% For every 1000 data points
for index = 1:1000:length(time)-1000
    idx1 = index;
    idx2 = index + 1000;
    % Compute the average frequency 
    [f(idx3),t(idx3)] = computeFrequency(time(idx1:idx2),accCh0(idx1:idx2));
    idx3 = idx3 + 1;
end
% Linear portion
f = f(35:775);
t = t(35:775);
% Linear fit
[coefficients, error] = polyfit(t,f,1);
[frequency, uncertainty] = polyval(coefficients,time,error);
% Frequency for fourier transform
f = Fs.*(0:L/2)/L;

%% FFT Tail
accCh1fft = fft(accCh1);
P2ch1 = abs(accCh1fft/L);
P1ch1 = P2ch1(1:L/2+1);
P1ch1(2:end-1) = 2*P1ch1(2:end-1);
P1ch1 = P1ch1/max(P1ch1);

%% FFT Wing
accCh2fft = fft(accCh2);
P2ch2 = abs(accCh2fft/L);
P1ch2 = P2ch2(1:L/2+1);
P1ch2(2:end-1) = 2*P1ch2(2:end-1);
P1ch2 = P1ch2/max(P1ch2);

%% FFT Nose
accCh3fft = fft(accCh3);
P2ch3 = abs(accCh3fft/L);
P1ch3 = P2ch3(1:L/2+1);
P1ch3(2:end-1) = 2*P1ch3(2:end-1);
P1ch3 = P1ch3/max(P1ch3);

%% FFT Vibrometer
accVibrofft = fft(vibro);
P2vibro = abs(accVibrofft/L);
P1vibro = P2vibro(1:L/2+1);
P1vibro(2:end-1) = 2*P1vibro(2:end-1);
P1vibro = P1vibro/max(P1vibro);

%% Excitation Frequency vs Time Plot
figure; hold on; grid minor;

plot(t,f);
plot(time, frequency, 'LineWidth',2)
plot(time, frequency+uncertainty,'k--')
plot(time, frequency-uncertainty,'k--')

xlabel("Time [s]");
ylabel("Frequency [Hz]")
legend("Frequency From Data", "Linear Regression")
title ("Excitation Frequency vs Time")

%% FFT Amplitude vs Excitation Frequency Plot 
figure; hold on; grid minor;

plot(f, P1ch1)
plot(f, P1ch2)
plot(f, P1ch3)
plot(f(636:end), P1vibro(636:end))

xlim([0 50]);
xlabel("Excitation Frequency [Hz]")
ylabel("FFT Amplitude")
legend("Tail", "Nose", "Wing", "Vibrometer")
title("Resonant Frequencies of Experimental Data")
