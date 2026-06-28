"""
Frequency Hopping Spread Spectrum (FHSS)

Spread spectrum is a technique for making a signal harder to intercept or jam
by spreading its energy across a wide band of frequencies. There are two main
flavors:

  - DSSS (Direct Sequence): spreads the signal across ALL frequencies at once
    by multiplying with a fast pseudo-random code.
  - FHSS (Frequency Hopping): uses only ONE narrow frequency at a time, but
    rapidly *hops* to a new frequency according to a pseudo-random pattern.

Think of FHSS like a conversation that keeps switching radio channels every
fraction of a second. An eavesdropper tuned to one channel only catches a
brief snippet before the signal jumps away. A jammer would need to blast
energy on EVERY possible channel simultaneously to be effective.

Key concepts:
  - Hopping pattern: the pseudo-random sequence that tells both transmitter
    and receiver which frequency to use next. Without it, you can't follow.
  - Hop rate: how many frequency changes per second.
  - Number of channels: how many distinct frequencies are available to hop
    among. More channels = harder to jam or intercept.

Real-world uses: Bluetooth (1600 hops/sec across 79 channels), military
radios, some Wi-Fi standards, and cordless phones.

Reference: converted from FHSS.m (Kashif Shahzad, 2004).
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.signal import spectrogram as _spectrogram


# ---------------------------------------------------------------------------
# Core functions (importable by notebook)
# ---------------------------------------------------------------------------

def generate_hopping_pattern(num_hops, num_channels, seed=None):
    """Generate a pseudo-random frequency hopping pattern.

    In real systems this pattern comes from a shared secret between
    transmitter and receiver (like a key). Here we just use a PRNG.

    Parameters
    ----------
    num_hops : int – how many hops (one per data symbol)
    num_channels : int – number of available frequency channels (1..num_channels)
    seed : int or None – random seed for reproducibility

    Returns
    -------
    pattern : ndarray of int – channel index for each hop (0-based)
    """
    rng = np.random.default_rng(seed)
    pattern = rng.integers(0, num_channels, size=num_hops)
    return pattern


def fhss_modulate(data, hop_freqs, fs, hop_rate, carrier_freq=0.0):
    """Apply FHSS modulation: BPSK on a carrier that hops frequencies.

    The process:
      1. Each data bit is BPSK-modulated (bit 0 -> -1, bit 1 -> +1).
      2. For each bit interval the carrier frequency is set according to the
         hopping pattern, so the signal jumps to a different frequency band.

    Parameters
    ----------
    data : array-like of 0/1 – the digital data bits
    hop_freqs : array-like of float – frequency (Hz) to use for each bit
        (length must equal len(data))
    fs : float – sample rate (samples/second)
    hop_rate : float – hops per second (determines samples per hop)
    carrier_freq : float – base carrier frequency offset (Hz). Each hop
        frequency is added on top of this.

    Returns
    -------
    t : ndarray – time vector for the full signal
    bpsk_signal : ndarray – BPSK modulated signal (before hopping)
    fhss_signal : ndarray – final frequency-hopped signal
    """
    data = np.asarray(data, dtype=int)
    hop_freqs = np.asarray(hop_freqs, dtype=float)
    num_bits = len(data)
    samples_per_hop = int(fs / hop_rate)
    total_samples = num_bits * samples_per_hop

    t = np.arange(total_samples) / fs

    # BPSK: map 0 -> -1, 1 -> +1
    bpsk_symbols = 2 * data - 1

    # Build the BPSK baseband signal (rectangular pulses)
    bpsk_signal = np.repeat(bpsk_symbols, samples_per_hop)

    # Build the hopping carrier: each hop interval uses a different frequency
    fhss_signal = np.zeros(total_samples)
    for i in range(num_bits):
        start = i * samples_per_hop
        end = start + samples_per_hop
        t_seg = t[start:end] - t[start]  # local time for phase continuity
        freq = carrier_freq + hop_freqs[i]
        fhss_signal[start:end] = bpsk_signal[start:end] * np.cos(
            2 * np.pi * freq * t_seg
        )

    return t, bpsk_signal, fhss_signal


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def plot_fhss_demo(num_bits=25, num_channels=6, hop_rate=100, fs=12000,
                   seed=42):
    """Generate the four-panel FHSS overview figure.

    Panels:
      1. BPSK modulated signal (data mapped to +1/-1)
      2. Hopping pattern (which channel is active over time — staircase)
      3. FHSS signal in the time domain
      4. Spectrogram showing the frequency hops (the key visual)
    """
    # --- Generate data and hopping pattern ---
    rng = np.random.default_rng(seed)
    data = rng.integers(0, 2, size=num_bits)

    pattern = generate_hopping_pattern(num_bits, num_channels, seed=seed)

    # Map channel indices to actual frequencies (evenly spaced)
    base_freq = 500       # lowest hop frequency (Hz)
    channel_spacing = 400  # Hz between channels
    channel_freqs = base_freq + np.arange(num_channels) * channel_spacing
    hop_freqs = channel_freqs[pattern]

    # --- Modulate ---
    t, bpsk_signal, fhss_signal = fhss_modulate(
        data, hop_freqs, fs, hop_rate
    )

    t_ms = t * 1000  # time in milliseconds

    # --- Build the staircase hopping pattern for plotting ---
    samples_per_hop = int(fs / hop_rate)
    hop_freq_trace = np.repeat(hop_freqs, samples_per_hop)

    # --- Plot ---
    fig, axes = plt.subplots(4, 1, figsize=(12, 10))
    fig.suptitle("Frequency Hopping Spread Spectrum (FHSS)", fontsize=14)

    # Panel 1: BPSK signal
    ax = axes[0]
    ax.plot(t_ms, bpsk_signal, color="tab:blue", linewidth=0.8)
    ax.set_title("BPSK Modulated Signal (bit 0 = -1, bit 1 = +1)")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude")
    ax.set_ylim(-1.5, 1.5)
    ax.grid(True, alpha=0.3)

    # Panel 2: Hopping pattern (staircase)
    ax = axes[1]
    ax.plot(t_ms, hop_freq_trace, color="tab:orange", linewidth=1.5,
            drawstyle="steps-post")
    ax.set_title("Hopping Pattern (pseudo-random channel sequence)")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Carrier Freq (Hz)")
    ax.set_yticks(channel_freqs)
    ax.grid(True, alpha=0.3)

    # Panel 3: FHSS time-domain signal
    ax = axes[2]
    ax.plot(t_ms, fhss_signal, color="tab:green", linewidth=0.4)
    ax.set_title("FHSS Signal (BPSK x hopping carrier)")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude")
    ax.set_ylim(-1.5, 1.5)
    ax.grid(True, alpha=0.3)

    # Panel 4: Spectrogram — the money plot
    ax = axes[3]
    nperseg = min(256, samples_per_hop)
    f, t_spec, Sxx = _spectrogram(fhss_signal, fs, nperseg=nperseg,
                                   noverlap=nperseg // 2)
    Sxx_db = 10 * np.log10(Sxx + 1e-12)
    im = ax.pcolormesh(t_spec * 1000, f, Sxx_db, shading="gouraud",
                       cmap="viridis")
    ax.set_title("Spectrogram — watch the signal hop between frequencies")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Frequency (Hz)")
    ax.set_ylim(0, channel_freqs[-1] + channel_spacing)
    fig.colorbar(im, ax=ax, label="Power (dB)")

    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Interactive widget (for Jupyter notebook)
# ---------------------------------------------------------------------------

def interactive_fhss():
    """Launch ipywidgets sliders to explore FHSS parameters.

    Sliders:
      - num_channels: how many frequency channels to hop among (2-12)
      - hop_rate: hops per second (50-500)
    """
    try:
        from ipywidgets import interact, IntSlider
    except ImportError:
        print("ipywidgets not available — run in Jupyter for interactive mode")
        return

    def _update(num_channels=6, hop_rate=100):
        fig = plot_fhss_demo(num_bits=25, num_channels=num_channels,
                             hop_rate=hop_rate, fs=12000, seed=42)
        plt.show()

    interact(
        _update,
        num_channels=IntSlider(value=6, min=2, max=12, step=1,
                               description="Channels"),
        hop_rate=IntSlider(value=100, min=50, max=500, step=25,
                           description="Hop Rate (Hz)"),
    )


# ---------------------------------------------------------------------------
# Standalone execution
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    matplotlib.use("Agg")
    fig = plot_fhss_demo()
    fig.savefig("fhss_demo.png", dpi=150)
    print("Saved fhss_demo.png")
