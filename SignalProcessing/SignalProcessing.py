import matplotlib.pyplot as plt
import numpy.random
from scipy import signal, fft

n = 500
Fs = 1000
F_max = 9
width = 21 / 2.54
height = 14 / 2.54
a = 0
b = 10
rnd = numpy.random.normal(a, b, n)
time_OX = numpy.arange(n)/Fs
w = F_max / (Fs / 2)
parameter = signal.butter(3, w, 'low', output='sos')
filtered_signal = signal.sosfiltfilt(parameter, rnd)
fig, ax = plt.subplots(figsize=(width, height))
ax.plot(time_OX, filtered_signal, linewidth="1")
ax.set_xlabel("Час (секунди)", fontsize=14)
ax.set_ylabel("Амплітуда сигналу", fontsize=14)
fig.savefig("./figures/Сигнал з максимальною частотою 29Гц.png")
spectrum = fft.fft(filtered_signal)
spectrum_abs = numpy.abs(fft.fftshift(spectrum))
frequency_counts = fft.fftfreq(500, 1 / 500)
frequency_counts_symetry = fft.fftshift(frequency_counts)
fig, ax = plt.subplots(figsize=(width, height))
ax.plot(frequency_counts_symetry, spectrum_abs, linewidth="1")
ax.set_xlabel("Частота", fontsize=14)
ax.set_ylabel("Амплітуда спектру", fontsize=14)
fig.savefig("./figures/Спектр сигналу з максимальною частотою 29Гц.png")
