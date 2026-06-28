"""
Signal Basics & FFT

Signals are just values that change over time. A radio wave, a sound wave, or
even a stock price are all signals. The most fundamental signal is a sinusoid
(sine wave), defined by three properties:
  - Frequency: how fast it oscillates (Hz = cycles per second)
  - Amplitude: how big the oscillations are
  - Phase: where in its cycle the wave starts

The key insight of DSP is that ANY signal can be broken down into a sum of
sinusoids at different frequencies. The FFT (Fast Fourier Transform) is the
algorithm that does this decomposition — it converts a signal from the "time
domain" (amplitude vs time) to the "frequency domain" (amplitude vs frequency).
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Core functions (importable by notebook)
# ---------------------------------------------------------------------------

def generate_sinusoid(freq, duration, fs, amplitude=1.0, phase=0.0):
    """Generate a sinusoidal signal.

    Parameters
    ----------
    freq : float – frequency in Hz
    duration : float – length in seconds
    fs : float – sample rate in samples/second
    amplitude : float – peak amplitude
    phase : float – phase offset in radians

    Returns
    -------
    t : ndarray – time vector
    signal : ndarray – sinusoidal samples
    """
    t = np.arange(0, duration, 1 / fs)
    signal = amplitude * np.sin(2 * np.pi * freq * t + phase)
    return t, signal


def compute_spectrum(signal, fs):
    """Compute single-sided amplitude spectrum via FFT.

    Returns
    -------
    freqs : ndarray – frequency bins (Hz)
    magnitude : ndarray – amplitude at each bin
    phase : ndarray – phase at each bin (radians)
    """
    N = len(signal)
    Y = np.fft.fft(signal)
    freqs = np.fft.fftfreq(N, 1 / fs)[:N // 2]
    magnitude = 2.0 / N * np.abs(Y[:N // 2])
    phase = np.angle(Y[:N // 2])
    return freqs, magnitude, phase


def compute_psd(signal, fs):
    """Compute power spectral density (PSD) in dB."""
    N = len(signal)
    Y = np.fft.fft(signal)
    psd = (1 / (fs * N)) * np.abs(Y[:N // 2]) ** 2
    psd[1:] *= 2
    freqs = np.fft.fftfreq(N, 1 / fs)[:N // 2]
    psd_db = 10 * np.log10(psd + 1e-12)
    return freqs, psd_db


def compute_spectrogram(signal, fs, nperseg=256):
    """Compute spectrogram using short-time FFT."""
    from scipy.signal import spectrogram as _spectrogram
    f, t, Sxx = _spectrogram(signal, fs, nperseg=nperseg)
    return f, t, 10 * np.log10(Sxx + 1e-12)


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def plot_signal_basics(fs=8000, duration=0.05):
    """Generate the four-panel overview figure."""
    f1, f2, f3 = 200, 500, 1300

    t, s1 = generate_sinusoid(f1, duration, fs, amplitude=1.0)
    _, s2 = generate_sinusoid(f2, duration, fs, amplitude=0.5)
    _, s3 = generate_sinusoid(f3, duration, fs, amplitude=0.3)
    composite = s1 + s2 + s3

    freqs, mag, phase = compute_spectrum(composite, fs)
    _, psd_db = compute_psd(composite, fs)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("Signal Basics & FFT", fontsize=14)

    # Time domain — individual sinusoids
    ax = axes[0, 0]
    ax.plot(t * 1000, s1, alpha=0.7, label=f"{f1} Hz")
    ax.plot(t * 1000, s2, alpha=0.7, label=f"{f2} Hz")
    ax.plot(t * 1000, s3, alpha=0.7, label=f"{f3} Hz")
    ax.set_title("Individual Sinusoids")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Time domain — composite
    ax = axes[0, 1]
    ax.plot(t * 1000, composite, color="tab:purple")
    ax.set_title("Composite Signal (sum of all three)")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude")
    ax.grid(True, alpha=0.3)

    # Frequency domain — magnitude spectrum
    ax = axes[1, 0]
    ax.stem(freqs, mag, linefmt="tab:blue", markerfmt="o", basefmt="k-")
    ax.set_title("FFT Magnitude Spectrum")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Amplitude")
    ax.set_xlim(0, 2000)
    ax.grid(True, alpha=0.3)

    # PSD
    ax = axes[1, 1]
    ax.plot(freqs, psd_db, color="tab:red")
    ax.set_title("Power Spectral Density")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Power (dB)")
    ax.set_xlim(0, 2000)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


def plot_spectrogram_demo(fs=8000, duration=0.5):
    """Show a spectrogram of a signal whose frequency changes over time."""
    t = np.arange(0, duration, 1 / fs)
    # Chirp: frequency sweeps from 200 Hz to 2000 Hz
    f0, f1 = 200, 2000
    signal = np.sin(2 * np.pi * (f0 * t + (f1 - f0) / (2 * duration) * t ** 2))

    f, t_spec, Sxx_db = compute_spectrogram(signal, fs, nperseg=128)

    fig, axes = plt.subplots(2, 1, figsize=(12, 6))
    fig.suptitle("Spectrogram — Frequency vs Time", fontsize=14)

    axes[0].plot(t * 1000, signal, color="tab:blue")
    axes[0].set_title("Chirp Signal (200 Hz -> 2000 Hz sweep)")
    axes[0].set_xlabel("Time (ms)")
    axes[0].set_ylabel("Amplitude")
    axes[0].grid(True, alpha=0.3)

    im = axes[1].pcolormesh(t_spec * 1000, f, Sxx_db, shading="gouraud", cmap="viridis")
    axes[1].set_title("Spectrogram")
    axes[1].set_xlabel("Time (ms)")
    axes[1].set_ylabel("Frequency (Hz)")
    fig.colorbar(im, ax=axes[1], label="Power (dB)")

    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Interactive widget (for Jupyter notebook)
# ---------------------------------------------------------------------------

def interactive_signal():
    """Launch ipywidgets sliders for exploring sinusoid superposition."""
    try:
        from ipywidgets import interact, FloatSlider
    except ImportError:
        print("ipywidgets not available — run in Jupyter for interactive mode")
        return

    def _update(f1=200, f2=500, f3=1300, a1=1.0, a2=0.5, a3=0.3):
        fs = 8000
        duration = 0.05
        t = np.arange(0, duration, 1 / fs)
        s = a1 * np.sin(2 * np.pi * f1 * t) + \
            a2 * np.sin(2 * np.pi * f2 * t) + \
            a3 * np.sin(2 * np.pi * f3 * t)

        freqs, mag, _ = compute_spectrum(s, fs)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
        ax1.plot(t * 1000, s)
        ax1.set_title("Time Domain")
        ax1.set_xlabel("Time (ms)")
        ax1.set_ylabel("Amplitude")
        ax1.grid(True, alpha=0.3)

        ax2.stem(freqs, mag, linefmt="tab:blue", markerfmt="o", basefmt="k-")
        ax2.set_title("Frequency Domain")
        ax2.set_xlabel("Frequency (Hz)")
        ax2.set_ylabel("Amplitude")
        ax2.set_xlim(0, 2000)
        ax2.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()

    freq_kw = dict(min=50, max=3000, step=50)
    amp_kw = dict(min=0.0, max=2.0, step=0.1)
    interact(
        _update,
        f1=FloatSlider(value=200, description="Freq 1 (Hz)", **freq_kw),
        f2=FloatSlider(value=500, description="Freq 2 (Hz)", **freq_kw),
        f3=FloatSlider(value=1300, description="Freq 3 (Hz)", **freq_kw),
        a1=FloatSlider(value=1.0, description="Amp 1", **amp_kw),
        a2=FloatSlider(value=0.5, description="Amp 2", **amp_kw),
        a3=FloatSlider(value=0.3, description="Amp 3", **amp_kw),
    )


# ---------------------------------------------------------------------------
# Standalone execution
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    matplotlib.use("Agg")
    fig1 = plot_signal_basics()
    fig1.savefig("signal_basics.png", dpi=150)
    print("Saved signal_basics.png")

    fig2 = plot_spectrogram_demo()
    fig2.savefig("spectrogram_demo.png", dpi=150)
    print("Saved spectrogram_demo.png")

    plt.show()
