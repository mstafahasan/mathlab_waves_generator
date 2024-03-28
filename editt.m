% Load the signal package
pkg load signal;

% Task 1: Generate signals for DO, RE, MI, FA
f0 = 440; % Fundamental frequency in Hz
alpha = 2^(1/12); % Alpha value for frequency calculation

% Define the frequencies for DO, RE, MI, and FA
f_DO = f0 * alpha^(-9);
f_RE = f0 * alpha^(-7);
f_MI = f0 * alpha^(-5);
f_FA = f0 * alpha^(-4);

% Time duration for each musical note (half a second)
duration = 0.5; % in seconds

% Generate time vector for each note
fs = 750; % Sampling frequency
% Time vector for each note
t_note = 0:1/fs:duration;

% Generate signals for each musical note
x1 = cos(2*pi*f_DO*t_note);
x2 = cos(2*pi*f_RE*t_note);
x3 = cos(2*pi*f_MI*t_note);
x4 = cos(2*pi*f_FA*t_note);

% Task 2: Create a signal sequentially playing DO, RE, MI, FA
x = [x1, x2, x3, x4];

filename = 'x(t)2nd.wav';

% Write the signal to a WAV file
audiowrite(filename, x, fs);

% Task 3: Plot the signal x(t) versus time t
t_total = 0:1/fs:(length(x)-1)/fs;
figure;
plot(t_total, x);
title('Signal x(t) versus Time');
xlabel('Time (seconds)');
ylabel('Amplitude');

% Task 4: Compute the energy of the signal x(t)
energy = sum(x.^2);
disp(['Energy of the signal x(t): ' num2str(energy)]);

% Task 5: Compute the frequency spectrum X(f) of the signal
N = length(x);
X = fft(x);
frequency_range = (-N/2:N/2-1) * (fs) / N; % Frequency range for FFT

% Task 6: Plot the magnitude of X(f) in the range -fs/2 <= f <= fs/2
figure;
plot(frequency_range, abs(fftshift(X)));
xlim([-2*fs, 2*fs]);
title('Magnitude of X(f)');
xlabel('Frequency (Hz)');
ylabel('|X(f)|');

% Task 7: Compute the Energy of the signal x(t) from its frequency spectrum X(f)
energy_from_spectrum = sum(abs(X).^2) / N;
disp(['Energy of the signal x(t) from its frequency spectrum: ' num2str(energy_from_spectrum)]);

% Task 8: Design a Butterworth low-pass filter
cutoff_frequency = f_RE; % Adjust as needed
filter_order = 20;
[b, a] = butter(filter_order, cutoff_frequency / (fs / 2), 'low');

% Task 9: Plot the magnitude and phase response of the Butterworth LPF
figure;
freqz(b, a);

% Task 10: Apply the signal x(t) to the Butterworth LPF
y1 = filter(b, a, x);

% Task 11: Store the generated signal y1(t) as a .wav file
filename_y1 = 'y1(t).wav';
audiowrite(filename_y1, y1, fs);

% Task 12: Plot the signal y1(t) versus time t
t_y1 = 0:1/fs:(length(y1)-1)/fs;
figure;
plot(t_y1, y1);
title('Signal y1(t) versus Time');
xlabel('Time (seconds)');
ylabel('Amplitude');

% Task 13: Compute the energy of the signal y1(t)
energy_y1 = sum(y1.^2);
disp(['Energy of the signal y1(t): ' num2str(energy_y1)]);

% Task 14: Compute the frequency spectrum Y1(f) of the signal
Y1 = fft(y1);
frequency_range_y1 = (-N/2:N/2-1) * (fs) / N; % Frequency range for FFT

% Task 15: Plot the magnitude of Y1(f) in the range -fs/2 <= f <= fs/2
figure;
plot(frequency_range_y1, abs(fftshift(Y1)));
xlim([-2*fs, 2*fs]);
title('Magnitude of Y1(f)');
xlabel('Frequency (Hz)');
ylabel('|Y1(f)|');

% Task 16: Compute the Energy of the signal y1(t) from its frequency spectrum Y1(f)
energy_from_spectrum_y1 = sum(abs(Y1).^2) / N;
disp(['Energy of the signal y1(t) from its frequency spectrum: ' num2str(energy_from_spectrum_y1)]);

% Task 17: Design a Butterworth high-pass filter
cutoff_frequency_hp = f_MI; % Adjust as needed
filter_order_hp = 20;
[b_hp, a_hp] = butter(filter_order_hp, cutoff_frequency_hp / (fs / 2), 'high');

% Task 18: Plot the magnitude and phase response of the Butterworth HPF
figure;
freqz(b_hp, a_hp);

% Task 19: Apply the signal x(t) to the Butterworth HPF
y2 = filter(b_hp, a_hp, x);

% Task 20: Store the generated signal y2(t) as a .wav file
filename_y2 = 'y2(t).wav';
audiowrite(filename_y2, y2, fs);

% Task 21: Plot the signal y2(t) versus time t
t_y2 = 0:1/fs:(length(y2)-1)/fs;
figure;
plot(t_y2, y2);
title('Signal y2(t) versus Time');
xlabel('Time (seconds)');
ylabel('Amplitude');

% Task 22: Compute the energy of the signal y2(t)
energy_y2 = sum(y2.^2);
disp(['Energy of the signal y2(t): ' num2str(energy_y2)]);

% Task 23: Compute the frequency spectrum Y2(f) of the signal
Y2 = fft(y2);
frequency_range_y2 = (-N/2:N/2-1) * (fs) / N; % Frequency range for FFT

% Task 24: Plot the magnitude of Y2(f) in the range -fs/2 <= f <= fs/2
figure;
plot(frequency_range_y2, abs(fftshift(Y2)));
xlim([-2*fs, 2*fs]);
title('Magnitude of Y2(f)');
xlabel('Frequency (Hz)');
ylabel('|Y2(f)|');

% Task 25: Compute the Energy of the signal y2(t) from its frequency spectrum Y2(f)
energy_from_spectrum_y2 = sum(abs(Y2).^2) / N;
disp(['Energy of the signal y2(t) from its frequency spectrum: ' num2str(energy_from_spectrum_y2)]);

