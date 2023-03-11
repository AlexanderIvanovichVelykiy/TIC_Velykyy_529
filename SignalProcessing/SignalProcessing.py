import numpy.random
from scipy import signal, fft
import matplotlib.pyplot as plt

def sg(a, b, n):
    return numpy.random.normal(a, b, n)
def determination(n, Fs):
    return numpy.arange(n) / Fs
def calculation_of_filter_parameters(F_max, Fs):
    w = F_max / (Fs / 2)
    return signal.butter(3, w, 'low', output='sos')
def signal_filtering(parameters_filter, sign):
    return signal.sosfiltfilt(parameters_filter, sign)
def calculation_of_the_signal_spectrum(f_sign, n):
    spectrum = fft.fft(f_sign)
    center_sp = numpy.abs(fft.fftshift(spectrum))
    return center_sp, ftt_reading(n)
def ftt_reading(n):
    f_readings = fft.fftfreq(n, 1 / n)
    return fft.fftshift(f_readings)
def tmpss(n, sign, F_filter, Fs):
    params = calculation_of_filter_parameters(F_filter, Fs)
    discrete_signals = []
    discrete_spectrums = []
    analog_signals = []
    dispersions = []
    signals_noise = []

    for Dt in [2, 4, 8, 16]:
        discrete_signal = numpy.zeros(n)
        for i in range(0, round(n / Dt)):
            discrete_signal[i * Dt] = sign[i * Dt]
        spectrum = fft.fft(discrete_signal)
        discrete_spectrums += [list(numpy.abs(fft.fftshift(spectrum)))]
        discrete_signals += [list(discrete_signal)]
        qwe = signal.sosfiltfilt(params, discrete_signal)
        e1 = qwe - sign
        dispersion = numpy.var(e1)
        dispersions += [dispersion]
        signals_noise += [numpy.var(sign)/dispersion]
        analog_signals += [list(qwe)]
    return discrete_signals, discrete_spectrums, analog_signals, dispersions, signals_noise

def dts(x, y, ox_txt, oy_txt, name):
    fig, ax = plt.subplots(figsize=(21 / 2.54, 14 / 2.54))
    ax.plot(x, y, linewidth=1)
    ax.set_xlabel(ox_txt, fontsize=14)
    ax.set_ylabel(oy_txt, fontsize=14)
    plt.title(name, fontsize=14)
    plt.grid()
    fig.savefig(f"./figures/{name}.png", dpi=600)

def ssd(x, y, ox_txt, oy_txt, name):
    fig, ax = plt.subplots(2, 2, figsize=(21 / 2.54, 14 / 2.54))
    s = 0
    for i in range(0, 2):
        for j in range(0, 2):
            ax[i][j].plot(x, y[s], linewidth=1)
            s += 1
    fig.supxlabel(ox_txt, fontsize=14)
    fig.supylabel(oy_txt, fontsize=14)
    fig.suptitle(name, fontsize=14)
    fig.savefig(f"./figures/{name}.png", dpi=600)

def sdts(x, y, ox_txt, oy_txt, name):
    fig, ax = plt.subplots(figsize=(21 / 2.54, 14 / 2.54))
    ax.plot(x, y, linewidth=1)
    ax.set_xlabel(ox_txt, fontsize=14)
    ax.set_ylabel(oy_txt, fontsize=14)
    plt.title(name, fontsize=14)
    plt.grid()
    fig.savefig(f"./figures/{name}.png", dpi=600)

def sssd(x, y, ox_txt, oy_txt, name):
        fig, ax = plt.subplots(2, 2, figsize=(21 / 2.54, 14 / 2.54))
        s = 0
        for i in range(0, 2):
            for j in range(0, 2):
                ax[i][j].plot(x, y[s], linewidth=1)
                s += 1
        fig.supxlabel(ox_txt, fontsize=14)
        fig.supylabel(oy_txt, fontsize=14)
        fig.suptitle(name, fontsize=14)
        fig.savefig(f"./figures/{name}.png", dpi=600)

def asd(x, y, ox_txt, oy_txt, name):
    fig, ax = plt.subplots(2, 2, figsize=(21 / 2.54, 14 / 2.54))
    s = 0
    for i in range(0, 2):
        for j in range(0, 2):
            ax[i][j].plot(x, y[s], linewidth=1)
            s += 1
    fig.supxlabel(ox_txt, fontsize=14)
    fig.supylabel(oy_txt, fontsize=14)
    fig.suptitle(name, fontsize=14)
    fig.savefig(f"./figures/{name}.png", dpi=600)

if __name__ == "__main__":
    n = 500
    Fs = 1000
    F_max = 9
    F_filter = 16
    sign = sg(0, 10, n)
    time_counts = determination(n, Fs)
    filter_params = calculation_of_filter_parameters(F_max, Fs)
    f_sign = signal_filtering(filter_params, sign)
    dis_sign, sp_sign, an_sign, dispersions, signals_noise = tmpss(n, f_sign, F_filter, Fs)
    xt = ftt_reading(n)
    ssd(time_counts, dis_sign, "Час (секунди)", "Амплітуда сигналу", "Сигнал з кроком дискретизації Dt = (2, 4, 8, 16)")
    sssd(xt, sp_sign, "Частота (Гц)", "Амплітуда спектру", "Спектри сигналів з кроком дискретизації Dt = (2, 4, 8, 16)")
    asd(time_counts, an_sign, "Час (секунди)", "Амплітуда сигналу", "Відновлені аналогові сигнали з кроком дискретизації Dt = (2, 4, 8, 16)")
    dts([2, 4, 8, 16], dispersions, "Крок дискретизації", "Дисперсія", "Залежність дисперсії від кроку дискретизації")
    sdts([2, 4, 8, 16], signals_noise, "Крок дискретизації", "ССШ", "Залежність співвідношення сигнал-шум від кроку дискретизації")
