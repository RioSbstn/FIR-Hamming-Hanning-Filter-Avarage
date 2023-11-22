% Path ke file Excel
file_path = '/MATLAB Drive/gyroscope.xls';

% Baca datasheet menggunakan readtable
datasheet = readtable(file_path);

% Akses data dari datasheet
t = datasheet{:, 1};
x = datasheet{:, 2};
y = datasheet{:, 3};
z = datasheet{:, 4};

% Plot data
figure(1)
subplot(3, 1, 1)
plot(t, x)
xlabel('Waktu')
ylabel('Amplitudo')
title('Sinyal X')

subplot(3, 1, 2)
plot(t, y)
xlabel('Waktu')
ylabel('Amplitudo')
title('Sinyal Y')

subplot(3, 1, 3)
plot(t, z)
xlabel('Waktu')
ylabel('Amplitudo')
title('Sinyal Z')

% Lakukan FFT pada sinyal
Fs = 1 / (t(2) - t(1));  % Frekuensi sampling (asumsi data diambil secara seragam)
X_fft = fft(x);
Y_fft = fft(y);
Z_fft = fft(z);
frekuensi = linspace(0, Fs, length(t));

% Plot hasil FFT
figure(2)
subplot(3, 1, 1)
plot(frekuensi, abs(X_fft))
xlabel('Frekuensi (Hz)')
ylabel('Magnitude')
title('FFT Sinyal X')

subplot(3, 1, 2)
plot(frekuensi, abs(Y_fft))
xlabel('Frekuensi (Hz)')
ylabel('Magnitude')
title('FFT Sinyal Y')

subplot(3, 1, 3)
plot(frekuensi, abs(Z_fft))
xlabel('Frekuensi (Hz)')
ylabel('Magnitude')
title('FFT Sinyal Z')

Fs = 500;
f = 10;
n = 1/Fs:1/Fs:1;
x = sin(2*pi*f*n);

%plot bandpass
fs = 1000;
b = fir1 (30,0.2);
[h,f]=freqz(b,1,512);
figure(3);
plot(f*fs/(2*pi),20*log10(abs(h)))
xlabel('Frekuensi');
ylabel('X axis');
title ('The gain response of bandpass filter')

%plot masukin filter
f1=400;
t=(0:200)/fs;
t1=(0.018:0.00001:12.769);
s=sin(2*f1*pi*t);
s1=sin(2*f1*pi*t1);
sf=filter(b,1,s);

figure (4)
subplot(2,1,1)
plot(t1,s1)
xlabel('time');
ylabel('x axis');
title('time domain diagram before filtering');


subplot(2,1,2)
Fs=fft(s,512);
AFs=abs(Fs);f=fs/512*(0:255);
plot(f,AFs(1:256));
xlabel('Frekuensi');
ylabel('X Axis');
title('Frekuensi domain sebelum Filtering');

figure(5)
subplot(2,1,1)
plot(t,sf)
xlabel('time');
ylabel('x axis');
title('time domain diagram after filtering');


subplot(2,1,2)
Fsf=fft(sf, 512);
AFsf=abs(Fsf);
f=(0:255)*fs/512;
plot(f,AFsf(1:256))
xlabel('Frequency');
ylabel('X Axis');
title('Frequency domain setelah filtering')

fs = 1000;
fc = 100;
m = 50;
N = 1024;
n = 0:1:m;

hdHn = fc/(fs) * sinc(fc/(fs) * (n - m/2));

wnHann = 0.5 * (1 - cos(2*pi*n/m));
hnHn = hdHn .* wnHann;
Hannfiltered = conv(x, hnHn);
fftFiltered = fft(Hannfiltered, N);

figure(6)
plot(Hannfiltered);
title('hann Filter');

hdHm = 0.2 * sinc(0.2 * (n - 6));
Hammfiltered = hdHm .* (0.54 - 0.46 * cos(2*pi*n/98));
fftHamm = fft(Hammfiltered, N);


figure(7)
plot(Hammfiltered);
title('Hamming Filter');


order = 2;
low_freq = 8.5;
high_freq = 11.5; % bandpass butterworth, passband antara 8.5 Hz dan 11.5 Hz.
Fs = 1000; % Frekuensi sampling yang sesuai
normalized_low_freq = low_freq / (Fs/2);
normalized_high_freq = high_freq / (Fs/2);
[b, a] = butter(order, [normalized_low_freq, normalized_high_freq], 'bandpass');

% Plot the magnitude response of the filter
[h, w] = freqz(b, a, 1024, Fs);
figure(8)
plot(w, 20*log10(abs(h))); % menghitung dan plotting magnitude respon dari tujuan filter
title('Magnitude Response of the Filter');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;



% Design Butterworth bandpass filter
order = 2;
low_freq = 8.5;
high_freq = 11.5; % bandpass butterworth, passband antara 8.5 Hz dan 11.5 Hz.
normalized_low_freq = low_freq / (Fs/2);
normalized_high_freq = high_freq / (Fs/2);
[b, a] = butter(order, [normalized_low_freq, normalized_high_freq], 'bandpass');

% Plot frequency response of the filter
figure(9)
freqz(b, a, 512, fs)
title('Frequency Response of Butterworth Bandpass Filter')

% spectral analysis filtered signal (frekuensi)
fs_filt = fft(y_filt, 512);
AFs_filt = abs(fs_filt);
f_filt = Fs/512 * (0:255);
figure(10)
subplot(2, 1, 2);
plot(f_filt, AFs_filt(1:256));
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Spectrum of Filtered Signal y(t)');

% Filter signal y
y_filt = filter(b, a, y); % Apply filter to the noise signal y, and store the result in y_filt

% Create a time vector for plotting
t_filt = linspace(0, (length(y_filt)-1)*(1/Fs), length(y_filt));

% Plot filtered signal
figure(11)
subplot(2, 1, 2)
plot(t_filt, y_filt)
title('Filtered Signal');
xlabel('Time (s)');
ylabel('Amplitude');

