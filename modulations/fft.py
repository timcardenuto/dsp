#!/usr/bin/env python

import numpy as np
import scipy as sp
from scipy import signal
from matplotlib import pyplot as plt
import sys


# Plotting
f = np.fromfile(open('/home/tc/Downloads/trimmedSamples.sigmf-data'), dtype=np.complex64)
x=f

# Frequency
X = 10*np.log10(np.abs(np.fft.fftshift(np.fft.fft(x)))**2)
f = np.linspace(-0.5, 0.5, len(X))
plt.subplot(2,2,1)
plt.plot(f, X)
plt.grid()
plt.xlabel("Frequency [Hz Normalized]")
plt.ylabel("PSD [dB]")

# Time
plt.subplot(2,2,2)
plt.plot(x.real[0:100])
plt.plot(x.imag[0:100])
plt.legend(['I','Q'])
plt.grid()
plt.xlabel("Time")
plt.ylabel("Sample")

# IQ
plt.subplot(2,2,4)
plt.plot(x.real[0:1000], x.imag[0:1000], '.')
plt.grid()
plt.xlabel("I")
plt.ylabel("Q")

# Spectogram
fft_size = 1024
sample_rate = 1e6
num_rows = int(np.floor(len(x)/fft_size))
spectrogram = np.zeros((num_rows, fft_size))
for i in range(num_rows):
    spectrogram[i,:] = 10*np.log10(np.abs(np.fft.fftshift(np.fft.fft(x[i*fft_size:(i+1)*fft_size])))**2)

plt.subplot(2,2,3)
plt.imshow(spectrogram, aspect='auto', extent = [sample_rate/-2, sample_rate/2, 0, len(x)/sample_rate])
plt.xlabel("Frequency [Hz]")
plt.ylabel("Time [s]")

plt.show()

sys.exit(0)


# f = np.fromfile(open('/home/tc/fm_101.8MHz_1Msps.cfile'), dtype=np.uint8)
f = np.fromfile(open('/home/tc/fm_101.8MHz_1Msps.cfile'), dtype=np.complex64)
center_freq = 101800000 # 101.8MHz
sample_rate = 1000000 # 1Msps
time_data = f[0:10000]

# create array of sample timestamps
num_samples = len(time_data) #len(f)
start_time = 0
end_time = num_samples / sample_rate  
timestamps = np.arange(start_time, end_time, (1/sample_rate))

# plot time series data
plt.subplot(2,2,1)
plt.plot(timestamps, time_data)
plt.xlabel('Sample Time')
plt.ylabel('Amplitude')


# get frequency data
freq_data = np.fft.fft(time_data)
baseband_freqs = np.fft.fftfreq(num_samples, (1/sample_rate))

#baseband_freqs_oneside = baseband_freqs[:num_samples//2]
#freq_data_oneside = np.abs(freq_data[:num_samples//2]) # i forget why we do this?

plt.subplot(2,2,2)
plt.plot(baseband_freqs, freq_data)
plt.title('Baseband Frequency Plot')
plt.xlabel('Freq (Hz)')
plt.ylabel('FFT Amplitude |X(freq)|')



# shift x-axis to RF, change units to MHz to make it easier to see
bandwidth = sample_rate / 2
half_bandwidth = bandwidth / 2
rf_freqs_oneside = baseband_freqs + (101800000-half_bandwidth)
rf_freqs_oneside_mhz = rf_freqs_oneside / 1000000

plt.subplot(2,2,4)
plt.plot(rf_freqs_oneside_mhz, freq_data)
plt.title('Radio Frequency Plot')
plt.xlabel('Freq (MHz)')
plt.ylabel('FFT Amplitude |X(freq)|')


f, t, Sxx = signal.spectrogram(time_data, sample_rate)
plt.subplot(2,2,3)
plt.pcolormesh(t, f, Sxx, shading='gouraud')
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [sec]')

plt.show()
