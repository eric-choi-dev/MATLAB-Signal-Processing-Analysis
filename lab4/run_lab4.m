%% Lab 4 Assignment Code 
clear;
clc;
close all; 

%% ========================================================================
%  Part A: The Fourier Transform and its Properties
%  ========================================================================

disp('Running Part A... (Opening all graphs)');


% Create the basic pulse x(t)
N = 100;
PulseWidth = 10;
t = 0:(N-1);
x = [ones(1, PulseWidth), zeros(1, N-PulseWidth)];

% Plot the waveform x(t)
figure;
stairs(t, x); 
grid on; 
axis([-10, 110, -0.1, 1.1]);
title('A: Original Pulse Signal x(t) [PulseWidth = 10]');
xlabel('Time (samples)');
ylabel('x(t)');

% Fourier Transform (FFT) of x(t)
Xf = fft(x);
f = [-(N/2):1:(N/2)-1]*(1/N); % Frequency vector

% Magnitude and Phase Spectrum of x(t)
figure;
subplot(2,1,1);
plot(f, fftshift(abs(Xf)));
grid on;
title('A: Magnitude Spectrum of x(t)');
ylabel('|X(f)|');
subplot(2,1,2);
plot(f, fftshift(angle(Xf)));
grid on;
title('A: Phase Spectrum of x(t)');
xlabel('Normalized Frequency');
ylabel('Angle(X(f))');

% --- A.1, A.2, A.4 ---
z_time = conv(x, x);
t_time = 0:(length(z_time)-1); 
N_conv = length(t_time); 
Xf_padded = fft(x, N_conv);
Zf_padded = Xf_padded .* Xf_padded; 
z_freq = ifft(Zf_padded);

% A.1 & A.4: Plot and compare z(t) results
figure;
subplot(2,1,1);
stairs(t_time, z_time);
grid on;
title('A.1/A.4: z(t) = x(t) * x(t) (Time-Domain conv)');
xlabel('Time (samples)');
ylabel('z(t)');
xlim([0, 30]);

subplot(2,1,2);
stairs(t_time, real(z_freq)); 
grid on;
title('A.4: z(t) (from Freq-Domain ifft(Xf .* Xf))');
xlabel('Time (samples)');
ylabel('z(t)');
xlim([0, 30]);

disp('A.4 Result: The two graphs should be identical (aside from numerical errors).');
disp('** Demonstrated Property: Convolution in the time domain is equal to multiplication in the frequency domain.**');

% --- A.3 ---
f_padded = [-(N_conv/2):1:(N_conv/2)-1]*(1/N_conv);
Zf_plot = fftshift(Zf_padded); 

figure;
subplot(2,1,1);
plot(f_padded, abs(Zf_plot));
grid on;
title('A.3: Magnitude Spectrum of Z(w)');
ylabel('|Z(f)|');
subplot(2,1,2);
plot(f_padded, angle(Zf_plot));
grid on;
title('A.3: Phase Spectrum of Z(w)');
xlabel('Normalized Frequency');
ylabel('Angle(Z(f))');


%% A.5: Time-Scaling Property
PulseWidth_5 = 5;
x_5 = [ones(1, PulseWidth_5), zeros(1, N-PulseWidth_5)];
Xf_5 = fft(x_5);
PulseWidth_25 = 25;
x_25 = [ones(1, PulseWidth_25), zeros(1, N-PulseWidth_25)];
Xf_25 = fft(x_25);

% Compare spectra
figure;
subplot(3,1,1);
plot(f, fftshift(abs(Xf))); grid on;
title('A.5: Magnitude Spectrum (Original PulseWidth = 10)');
ylabel('|X(f)|');
subplot(3,1,2);
plot(f, fftshift(abs(Xf_5))); grid on;
title('A.5: Magnitude Spectrum (PulseWidth = 5)');
ylabel('|X_5(f)|');
subplot(3,1,3);
plot(f, fftshift(abs(Xf_25))); grid on;
title('A.5: Magnitude Spectrum (PulseWidth = 25)');
xlabel('Normalized Frequency');
ylabel('|X_{25}(f)|');

disp('A.5 Result: When the pulse is narrower in time (Width=5), its spectrum is wider in frequency.');
disp('Conversely, when the pulse is wider in time (Width=25), its spectrum is narrower.');
disp('** Demonstrated Property: Time-Scaling Property.**');


%% A.6: Frequency-Shifting (Modulation) Property
t = 0:(N-1); 
w_plus = x .* exp(1j * (pi/3) .* t);   
w_minus = x .* exp(-1j * (pi/3) .* t); 
w_c = x .* cos((pi/3) .* t);        
W_plus_f = fft(w_plus);
W_minus_f = fft(w_minus);
W_c_f = fft(w_c);

% Plot spectra
figure;
subplot(3,1,1);
plot(f, fftshift(abs(W_plus_f))); grid on;
title('A.6: Mag Spectrum of w_+(t) = x(t) * e^{j(\pi/3)t}');
ylabel('|W_+(f)|');
subplot(3,1,2);
plot(f, fftshift(abs(W_minus_f))); grid on;
title('A.6: Mag Spectrum of w_-(t) = x(t) * e^{-j(\pi/3)t}');
ylabel('|W_-(f)|');
subplot(3,1,3);
plot(f, fftshift(abs(W_c_f))); grid on;
title('A.6: Mag Spectrum of w_c(t) = x(t) * cos(\pi/3)t');
xlabel('Normalized Frequency');
ylabel('|W_c(f)|');

disp('A.6 Result: The original spectrum (A.3) is shifted in frequency by multiplication with exp(j...).');
disp('Multiplication by cos(...) shifts the spectrum to both positive and negative frequencies.');
disp('** Demonstrated Property: Frequency-Shifting (Modulation) Property.**');


%% ========================================================================
%  Part B: Application of the Fourier Transform (Modulation System)
%  ========================================================================

disp('Running Part B... (Opening all graphs)');

% --- B.1 ---

% 1. Load data and set Fs
load('Lab4_Data.mat'); 
Fs = 32000; 
N = length(xspeech); 

% 2. Analyze Signal and Channel Characteristics
figure;
subplot(2,2,1);
MagSpect(xspeech);
title('B.1: Original Speech Signal Spectrum (xspeech)');
subplot(2,2,2);
MagSpect(hLPF2000);
title('B.1: LPF (2.0 kHz) Spectrum');
subplot(2,2,3);
MagSpect(hLPF2500);
title('B.1: LPF (2.5 kHz) Spectrum');
subplot(2,2,4);
MagSpect(hChannel);
title('B.1: Transmission Channel Spectrum (hChannel)');

disp('B.1: Analyzed signal and channel spectra.');
disp('The speech signal (0-3.5kHz) cannot pass through the channel (8-12kHz).');
disp('Solution: Filter and Modulate the speech signal to shift it into the channel band.');

% 3. Transmitter (Coder) Design
disp('Designing Transmitter (Coder)...');
x_filtered = conv(xspeech, hLPF2000, 'same');
Fc = 6500; % 6.5kHz carrier frequency
carrier = osc(Fc, N, Fs); 
s_t = x_filtered .* carrier;

% 4. Pass through Channel
disp('Transmitting signal through channel...');
r_t = conv(s_t, hChannel, 'same');

% 5. Receiver (Decoder) Design
disp('Designing Receiver (Decoder)...');
v_t = r_t .* carrier;
x_recovered = conv(v_t, hLPF2500, 'same');

% 6. Compare Resulting Spectra
figure;
subplot(2,2,1);
MagSpect(x_filtered);
title('B.1: Filtered Speech (before modulation)');
subplot(2,2,2);
MagSpect(s_t);
title('B.1: Modulated Signal (before channel)');
subplot(2,2,3);
MagSpect(v_t);
title('B.1: Demodulated Signal (before LPF)');
subplot(2,2,4);
MagSpect(x_recovered);
title('B.1: Final Recovered Speech');

disp('B.1: Design complete. Check the spectra.');

% 7. Listen to Audio Signals (Comparison)
disp('Playing audio signals... (Check volume)');

disp('Playing: Original Speech (xspeech)');
sound(xspeech, Fs);
pause(3); 

disp('Playing: Filtered Speech (x_filtered) - 2kHz LPF');
sound(x_filtered, Fs);
pause(3);

disp('Playing: Final Recovered Speech (x_recovered)');
sound(x_recovered / max(abs(x_recovered)), Fs);
pause(3);

disp('Lab 4 (Parts A & B) is complete.');