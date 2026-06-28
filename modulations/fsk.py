"""
Frequency Shift Keying (FSK)

FSK encodes digital data by switching between different carrier frequencies.
Think of it like Morse code, but instead of dots and dashes you use two
(or more) musical notes: one note means "0", a different note means "1".

Binary FSK (the simplest form) uses exactly two frequencies:
  - f0: the frequency transmitted when the data bit is 0
  - f1: the frequency transmitted when the data bit is 1

Each bit is held for a fixed "symbol duration", so the transmitter hops
between f0 and f1 as it reads through the bit sequence.

Key concepts:
  - Tone spacing: the gap between f0 and f1 (|f1 - f0|). Wider spacing
    makes the two tones easier to tell apart at the receiver, but uses
    more bandwidth.
  - Orthogonality: if the tone spacing equals 1/(2*symbol_duration), the
    two tones are mathematically orthogonal — their correlation is zero,
    which gives the best detection performance for a given spacing.
  - Comparison to ASK/PSK: ASK varies amplitude, PSK varies phase, and
    FSK varies frequency. FSK wastes more bandwidth than PSK (less
    spectrally efficient), but because amplitude stays constant, FSK is
    very robust to amplitude fading — useful in noisy channels like
    HF radio, pagers, and low-power IoT (e.g. LoRa uses a chirped FSK).
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.signal import spectrogram as _spectrogram


# ---------------------------------------------------------------------------
# Core functions (importable by notebooks)
# ---------------------------------------------------------------------------

def fsk_modulate(bits, f0, f1, fs, symbol_duration):
    """Generate a binary-FSK waveform from a sequence of bits.

    For each bit the function produces a sinusoid at f0 (bit=0) or f1
    (bit=1) lasting *symbol_duration* seconds, then concatenates them.

    Parameters
    ----------
    bits : array-like – sequence of 0s and 1s
    f0 : float – carrier frequency for bit 0 (Hz)
    f1 : float – carrier frequency for bit 1 (Hz)
    fs : float – sample rate (samples/second)
    symbol_duration : float – duration of each bit in seconds

    Returns
    -------
    t : ndarray – time vector for the full waveform
    signal : ndarray – FSK-modulated samples
    """
    bits = np.asarray(bits, dtype=int)
    samples_per_symbol = int(fs * symbol_duration)
    t_symbol = np.arange(samples_per_symbol) / fs

    # Pre-compute one period of each tone (for a single symbol duration)
    tone0 = np.sin(2 * np.pi * f0 * t_symbol)
    tone1 = np.sin(2 * np.pi * f1 * t_symbol)

    # Build the modulated signal by selecting the appropriate tone per bit
    signal = np.concatenate([tone1 if b else tone0 for b in bits])
    t = np.arange(len(signal)) / fs
    return t, signal


def fsk_demodulate(signal, f0, f1, fs, symbol_duration):
    """Demodulate a binary-FSK waveform back into bits.

    Uses energy detection: for each symbol window, correlate the received
    signal with reference sinusoids at f0 and f1. Whichever correlation
    has more energy wins — that determines the bit value.

    Parameters
    ----------
    signal : ndarray – received FSK waveform
    f0 : float – frequency representing bit 0 (Hz)
    f1 : float – frequency representing bit 1 (Hz)
    fs : float – sample rate (samples/second)
    symbol_duration : float – duration of each bit in seconds

    Returns
    -------
    bits : ndarray – recovered bit sequence (0s and 1s)
    """
    samples_per_symbol = int(fs * symbol_duration)
    t_symbol = np.arange(samples_per_symbol) / fs
    ref0 = np.sin(2 * np.pi * f0 * t_symbol)
    ref1 = np.sin(2 * np.pi * f1 * t_symbol)

    num_symbols = len(signal) // samples_per_symbol
    bits = np.zeros(num_symbols, dtype=int)

    for i in range(num_symbols):
        start = i * samples_per_symbol
        end = start + samples_per_symbol
        chunk = signal[start:end]

        # Correlate with each reference tone — higher energy wins
        energy0 = np.abs(np.sum(chunk * ref0))
        energy1 = np.abs(np.sum(chunk * ref1))
        bits[i] = 1 if energy1 > energy0 else 0

    return bits


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def plot_fsk_demo(bits=None, f0=200, f1=600, fs=8000, symbol_duration=0.01):
    """Generate a multi-panel figure illustrating FSK modulation.

    Panels
    ------
    1. Original bit sequence (stem plot)
    2. FSK modulated waveform — you can see the frequency change per bit
    3. Spectrogram — a heat-map of frequency vs time showing the hops
    4. Individual tone spectra — magnitude spectrum of one symbol at f0
       and one at f1, proving they sit at distinct frequencies
    """
    if bits is None:
        bits = [1, 0, 1, 1, 0, 0, 1, 0]
    bits = np.asarray(bits, dtype=int)

    t, signal = fsk_modulate(bits, f0, f1, fs, symbol_duration)

    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    fig.suptitle("Frequency Shift Keying (FSK)", fontsize=14)

    # --- Panel 1: original bit sequence ---
    ax = axes[0, 0]
    bit_times = np.arange(len(bits)) * symbol_duration * 1000  # ms
    ax.step(bit_times, bits, where="post", color="tab:blue", linewidth=2)
    for i, b in enumerate(bits):
        ax.annotate(str(b), (bit_times[i] + symbol_duration * 500, b),
                    ha="center", va="bottom", fontsize=9, color="tab:blue")
    ax.set_ylim(-0.3, 1.5)
    ax.set_title("Original Bit Sequence")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Bit Value")
    ax.grid(True, alpha=0.3)

    # --- Panel 2: FSK modulated signal ---
    ax = axes[0, 1]
    ax.plot(t * 1000, signal, color="tab:orange", linewidth=0.6)
    # Draw vertical lines at bit boundaries
    for i in range(1, len(bits)):
        ax.axvline(i * symbol_duration * 1000, color="gray",
                   linestyle="--", alpha=0.4)
    ax.set_title(f"FSK Modulated Signal  (f0={f0} Hz, f1={f1} Hz)")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude")
    ax.grid(True, alpha=0.3)

    # --- Panel 3: spectrogram showing frequency hops ---
    ax = axes[1, 0]
    # Choose nperseg so we get reasonable time/frequency resolution
    nperseg = min(len(signal), int(fs * symbol_duration) // 2)
    nperseg = max(nperseg, 32)  # floor to avoid degenerate FFTs
    f_spec, t_spec, Sxx = _spectrogram(signal, fs, nperseg=nperseg,
                                       noverlap=nperseg * 3 // 4)
    Sxx_db = 10 * np.log10(Sxx + 1e-12)
    im = ax.pcolormesh(t_spec * 1000, f_spec, Sxx_db,
                       shading="gouraud", cmap="viridis")
    ax.set_ylim(0, max(f0, f1) * 2.5)
    ax.set_title("Spectrogram — Frequency Hops Over Time")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Frequency (Hz)")
    fig.colorbar(im, ax=ax, label="Power (dB)")

    # --- Panel 4: individual tone spectra ---
    ax = axes[1, 1]
    samples_per_symbol = int(fs * symbol_duration)
    t_sym = np.arange(samples_per_symbol) / fs
    tone0 = np.sin(2 * np.pi * f0 * t_sym)
    tone1 = np.sin(2 * np.pi * f1 * t_sym)

    N = len(tone0)
    freqs = np.fft.fftfreq(N, 1 / fs)[:N // 2]
    mag0 = 2.0 / N * np.abs(np.fft.fft(tone0)[:N // 2])
    mag1 = 2.0 / N * np.abs(np.fft.fft(tone1)[:N // 2])

    ax.plot(freqs, mag0, label=f"Tone for bit 0 ({f0} Hz)",
            color="tab:green", linewidth=1.5)
    ax.plot(freqs, mag1, label=f"Tone for bit 1 ({f1} Hz)",
            color="tab:red", linewidth=1.5)
    ax.set_xlim(0, max(f0, f1) * 2.5)
    ax.set_title("Individual Tone Spectra")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Amplitude")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Interactive widget (for Jupyter notebook)
# ---------------------------------------------------------------------------

def interactive_fsk():
    """Launch ipywidgets sliders for exploring FSK modulation.

    Sliders let you adjust f0, f1, and the symbol rate, then instantly
    see how the modulated waveform and spectrogram change.
    """
    try:
        from ipywidgets import interact, FloatSlider, IntSlider
    except ImportError:
        print("ipywidgets not available — run in Jupyter for interactive mode")
        return

    bits = [1, 0, 1, 1, 0, 0, 1, 0]
    fs = 8000

    def _update(f0=200, f1=600, symbol_rate=100):
        symbol_duration = 1.0 / symbol_rate
        t, signal = fsk_modulate(bits, f0, f1, fs, symbol_duration)

        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 8))
        fig.suptitle("Interactive FSK Explorer", fontsize=13)

        # Bit sequence as digital waveform
        samples_per_sym = int(fs * symbol_duration)
        digital = np.concatenate([np.full(samples_per_sym, b) for b in bits])
        ax1.plot(t * 1000, digital, color="tab:blue", linewidth=1.5)
        ax1.set_title("Digital Bit Stream")
        ax1.set_ylabel("Bit Value")
        ax1.set_ylim(-0.3, 1.5)
        ax1.grid(True, alpha=0.3)

        # FSK waveform
        ax2.plot(t * 1000, signal, color="tab:orange", linewidth=0.6)
        for i in range(1, len(bits)):
            ax2.axvline(i * symbol_duration * 1000, color="gray",
                        linestyle="--", alpha=0.4)
        ax2.set_title(f"FSK Signal  (f0={f0} Hz, f1={f1} Hz, "
                      f"rate={symbol_rate} sym/s)")
        ax2.set_ylabel("Amplitude")
        ax2.grid(True, alpha=0.3)

        # Spectrogram
        nperseg = min(len(signal), samples_per_sym // 2)
        nperseg = max(nperseg, 32)
        f_spec, t_spec, Sxx = _spectrogram(signal, fs, nperseg=nperseg,
                                           noverlap=nperseg * 3 // 4)
        Sxx_db = 10 * np.log10(Sxx + 1e-12)
        ax3.pcolormesh(t_spec * 1000, f_spec, Sxx_db,
                       shading="gouraud", cmap="viridis")
        ax3.set_ylim(0, max(f0, f1) * 2.5)
        ax3.set_title("Spectrogram")
        ax3.set_xlabel("Time (ms)")
        ax3.set_ylabel("Frequency (Hz)")
        ax3.grid(False)

        plt.tight_layout()
        plt.show()

    interact(
        _update,
        f0=FloatSlider(value=200, min=50, max=2000, step=50,
                       description="f0 (Hz)"),
        f1=FloatSlider(value=600, min=50, max=2000, step=50,
                       description="f1 (Hz)"),
        symbol_rate=IntSlider(value=100, min=10, max=500, step=10,
                              description="Sym Rate"),
    )


# ---------------------------------------------------------------------------
# Standalone execution
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    matplotlib.use("Agg")

    fig = plot_fsk_demo()
    fig.savefig("fsk_demo.png", dpi=150)
    print("Saved fsk_demo.png")

    plt.show()
