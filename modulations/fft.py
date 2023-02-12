#!/usr/bin/env python

import numpy as np
import scipy as sp
from scipy import signal
from matplotlib import pyplot as plt


f = np.fromfile(open('/home/tc/fm_101.8MHz_1Msps.cfile'), dtype=np.uint8)
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


f, t, Sxx = signal.spectrogram(time_data, center_freq)
plt.subplot(2,2,3)
plt.pcolormesh(t, f, Sxx, shading='gouraud')
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [sec]')

plt.show()
plt.show()
