"""
BPSK — Binary Phase Shift Keying

BPSK is the simplest form of digital modulation. The idea is straightforward:
transmit data by flipping the phase of a carrier wave.

  - To send a binary 1, transmit the carrier as-is        (phase = 0 degrees).
  - To send a binary 0, flip the carrier upside down       (phase = 180 degrees).

The receiver just checks whether the incoming wave is right-side-up or upside
down to recover each bit. Because there are only two states, the "constellation
diagram" (a scatter plot of received signal samples on the I/Q plane) shows two
clusters: one at +1 on the real axis (bit 1) and one at -1 (bit 0).

Why is this useful?  BPSK is the most noise-resistant modulation — the two
constellation points are as far apart as possible for a given transmit power.
The trade-off is speed: each symbol only carries one bit.

Pulse shaping (RRC — Root Raised Cosine) is applied after mapping bits to
symbols. Without it, the abrupt transitions between symbols would spread the
signal across a wide bandwidth.  The RRC filter smooths these transitions so
the signal fits neatly into its allocated frequency band, without causing
interference to neighboring channels.

BER (Bit Error Rate) measures link quality: the fraction of received bits that
are wrong.  For BPSK in additive white Gaussian noise, BER = Q(sqrt(2*Eb/N0)),
where Eb/N0 is the energy-per-bit to noise ratio.
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Core functions (importable by notebook)
# ---------------------------------------------------------------------------

def bpsk_modulate(bits, carrier_freq, fs, samples_per_symbol):
    """Map a bit sequence to a BPSK waveform.

    Each bit is mapped to a symbol (+1 for bit 1, -1 for bit 0), held for
    *samples_per_symbol* samples, then multiplied by a cosine carrier.

    Parameters
    ----------
    bits : array-like of int – sequence of 0s and 1s
    carrier_freq : float – carrier frequency in Hz
    fs : float – sample rate in samples/second
    samples_per_symbol : int – how many samples represent one bit

    Returns
    -------
    t : ndarray – time vector (seconds)
    modulated : ndarray – BPSK waveform
    baseband : ndarray – NRZ baseband signal (+1 / -1) before carrier mixing
    """
    bits = np.asarray(bits, dtype=int)
    # Map bits to NRZ symbols: 1 -> +1, 0 -> -1
    symbols = 2 * bits - 1
    # Repeat each symbol for samples_per_symbol samples
    baseband = np.repeat(symbols, samples_per_symbol)
    num_samples = len(baseband)
    t = np.arange(num_samples) / fs
    carrier = np.cos(2 * np.pi * carrier_freq * t)
    modulated = baseband * carrier
    return t, modulated, baseband


def bpsk_demodulate(signal, carrier_freq, fs, samples_per_symbol):
    """Recover bits from a BPSK waveform using coherent demodulation.

    Multiplies the received signal by a local copy of the carrier (coherent
    detection), then integrates over each symbol period.  The sign of the
    integral decides whether the bit is 0 or 1.

    Parameters
    ----------
    signal : ndarray – received BPSK waveform
    carrier_freq : float – carrier frequency in Hz (must match transmitter)
    fs : float – sample rate in samples/second
    samples_per_symbol : int – samples per bit (must match transmitter)

    Returns
    -------
    bits : ndarray of int – recovered bit sequence
    """
    num_samples = len(signal)
    t = np.arange(num_samples) / fs
    carrier = np.cos(2 * np.pi * carrier_freq * t)
    # Multiply by carrier — this moves the baseband content back to DC
    mixed = signal * carrier
    # Integrate-and-dump: sum samples within each symbol period
    num_symbols = num_samples // samples_per_symbol
    bits = np.zeros(num_symbols, dtype=int)
    for i in range(num_symbols):
        chunk = mixed[i * samples_per_symbol:(i + 1) * samples_per_symbol]
        bits[i] = 1 if np.sum(chunk) > 0 else 0
    return bits


def rrc_filter(num_taps, sps, beta):
    """Design a Root Raised Cosine (RRC) pulse-shaping filter.

    The RRC filter is applied at the transmitter to shape each symbol's pulse
    so that the transmitted spectrum is compact.  When the same filter is
    applied at the receiver, the combined response is a Raised Cosine — which
    has zero inter-symbol interference at the ideal sampling instants.

    Parameters
    ----------
    num_taps : int – filter length (odd is typical, e.g. 101)
    sps : int – samples per symbol (oversampling factor)
    beta : float – roll-off factor, 0 to 1 (0 = brickwall, 1 = widest)

    Returns
    -------
    h : ndarray – filter coefficients, normalized to unit energy
    """
    t = np.arange(num_taps) - (num_taps - 1) // 2
    h = np.zeros(num_taps)
    for i, ti in enumerate(t):
        if ti == 0:
            h[i] = 1.0 / sps * (1 + beta * (4 / np.pi - 1))
        elif abs(ti) == sps / (4 * beta + 1e-30):
            h[i] = beta / (sps * np.sqrt(2)) * (
                (1 + 2 / np.pi) * np.sin(np.pi / (4 * beta))
                + (1 - 2 / np.pi) * np.cos(np.pi / (4 * beta))
            )
        else:
            num = np.sin(np.pi * ti / sps * (1 - beta)) + \
                  4 * beta * ti / sps * np.cos(np.pi * ti / sps * (1 + beta))
            den = np.pi * ti / sps * (1 - (4 * beta * ti / sps) ** 2)
            h[i] = num / (den + 1e-30)
    # Normalize to unit energy
    h /= np.sqrt(np.sum(h ** 2))
    return h


def add_awgn(signal, snr_db):
    """Add white Gaussian noise to achieve a target SNR (in dB).

    Parameters
    ----------
    signal : ndarray – clean signal
    snr_db : float – desired signal-to-noise ratio in decibels

    Returns
    -------
    noisy : ndarray – signal + noise
    """
    sig_power = np.mean(signal ** 2)
    noise_power = sig_power / (10 ** (snr_db / 10))
    noise = np.sqrt(noise_power) * np.random.randn(len(signal))
    return signal + noise


def compute_psd(signal, fs):
    """Compute power spectral density (PSD) in dB."""
    N = len(signal)
    Y = np.fft.fft(signal)
    psd = (1 / (fs * N)) * np.abs(Y[:N // 2]) ** 2
    psd[1:] *= 2
    freqs = np.fft.fftfreq(N, 1 / fs)[:N // 2]
    psd_db = 10 * np.log10(psd + 1e-12)
    return freqs, psd_db


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def plot_bpsk_demo(num_bits=12, carrier_freq=100, fs=4000,
                   samples_per_symbol=120, snr_db=20):
    """Generate a multi-panel figure illustrating BPSK modulation.

    Panels
    ------
    1. Original bit sequence (stem plot)
    2. BPSK modulated waveform with bit boundaries
    3. IQ constellation diagram (scatter plot of received symbols)
    4. Power spectral density of the modulated signal
    """
    np.random.seed(42)
    bits = np.random.randint(0, 2, num_bits)

    t, modulated, baseband = bpsk_modulate(bits, carrier_freq, fs,
                                           samples_per_symbol)
    noisy = add_awgn(modulated, snr_db)

    # Demodulate to get constellation points (correlator outputs per symbol)
    t_vec = np.arange(len(noisy)) / fs
    carrier_local = np.cos(2 * np.pi * carrier_freq * t_vec)
    mixed = noisy * carrier_local
    iq_i = np.array([np.mean(mixed[i * samples_per_symbol:(i + 1) * samples_per_symbol])
                     for i in range(num_bits)])
    # For pure BPSK the Q component is ideally zero; show noise on Q axis
    carrier_q = -np.sin(2 * np.pi * carrier_freq * t_vec)
    mixed_q = noisy * carrier_q
    iq_q = np.array([np.mean(mixed_q[i * samples_per_symbol:(i + 1) * samples_per_symbol])
                     for i in range(num_bits)])

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle("BPSK — Binary Phase Shift Keying", fontsize=14)

    # --- Panel 1: Original bits ---
    ax = axes[0, 0]
    bit_indices = np.arange(num_bits)
    ax.stem(bit_indices, bits, linefmt="tab:blue", markerfmt="o",
            basefmt="k-")
    ax.set_title("Original Bit Sequence")
    ax.set_xlabel("Bit index")
    ax.set_ylabel("Bit value")
    ax.set_ylim(-0.3, 1.5)
    ax.set_yticks([0, 1])
    ax.grid(True, alpha=0.3)

    # --- Panel 2: Modulated waveform ---
    ax = axes[0, 1]
    t_ms = t * 1000
    ax.plot(t_ms, modulated, color="tab:blue", linewidth=0.6)
    # Draw bit boundaries
    for i in range(1, num_bits):
        boundary = i * samples_per_symbol / fs * 1000
        ax.axvline(boundary, color="tab:red", linewidth=0.5, linestyle="--",
                   alpha=0.5)
    # Label each bit above the waveform
    for i, b in enumerate(bits):
        center = (i + 0.5) * samples_per_symbol / fs * 1000
        ax.text(center, 1.15, str(b), ha="center", fontsize=7,
                color="tab:red")
    ax.set_title("BPSK Modulated Waveform")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude")
    ax.set_ylim(-1.5, 1.5)
    ax.grid(True, alpha=0.3)

    # --- Panel 3: Constellation diagram ---
    ax = axes[1, 0]
    colors = ["tab:red" if b == 0 else "tab:blue" for b in bits]
    ax.scatter(iq_i, iq_q, c=colors, s=60, edgecolors="k", linewidths=0.5,
               zorder=3)
    ax.axhline(0, color="k", linewidth=0.5)
    ax.axvline(0, color="k", linewidth=0.5)
    ax.set_title("IQ Constellation Diagram")
    ax.set_xlabel("In-phase (I)")
    ax.set_ylabel("Quadrature (Q)")
    ax.set_aspect("equal")
    # Draw ideal constellation points
    ax.plot([-0.5, 0.5], [0, 0], "x", color="gray", markersize=12,
            markeredgewidth=2, label="Ideal points")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.3)

    # --- Panel 4: Power spectral density ---
    ax = axes[1, 1]
    freqs, psd_db = compute_psd(modulated, fs)
    ax.plot(freqs, psd_db, color="tab:purple")
    ax.set_title("Power Spectral Density")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Power (dB)")
    ax.set_xlim(0, fs / 2)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Interactive widget (for Jupyter notebook)
# ---------------------------------------------------------------------------

def interactive_bpsk():
    """Launch ipywidgets sliders to explore BPSK under varying noise.

    Sliders
    -------
    - SNR (dB): controls how much noise is added. Lower values smear the
      constellation points and increase bit errors.
    - Number of bits: controls how many symbols appear in the constellation.
    """
    try:
        from ipywidgets import interact, IntSlider, FloatSlider
    except ImportError:
        print("ipywidgets not available — run in Jupyter for interactive mode")
        return

    def _update(snr_db=20.0, num_bits=50):
        carrier_freq = 100
        fs = 4000
        sps = 120

        np.random.seed(0)
        bits = np.random.randint(0, 2, num_bits)
        t, modulated, baseband = bpsk_modulate(bits, carrier_freq, fs, sps)
        noisy = add_awgn(modulated, snr_db)

        # Recover bits
        recovered = bpsk_demodulate(noisy, carrier_freq, fs, sps)

        # Compute constellation points
        t_vec = np.arange(len(noisy)) / fs
        carrier_i = np.cos(2 * np.pi * carrier_freq * t_vec)
        carrier_q = -np.sin(2 * np.pi * carrier_freq * t_vec)
        mixed_i = noisy * carrier_i
        mixed_q = noisy * carrier_q
        iq_i = np.array([np.mean(mixed_i[k * sps:(k + 1) * sps])
                         for k in range(num_bits)])
        iq_q = np.array([np.mean(mixed_q[k * sps:(k + 1) * sps])
                         for k in range(num_bits)])

        # BER
        min_len = min(len(bits), len(recovered))
        errors = np.sum(bits[:min_len] != recovered[:min_len])
        ber = errors / min_len

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        fig.suptitle(f"BPSK Interactive  |  SNR = {snr_db:.1f} dB  |  "
                     f"BER = {ber:.4f}  ({errors}/{min_len} errors)",
                     fontsize=12)

        # Constellation
        colors = ["tab:red" if b == 0 else "tab:blue" for b in bits]
        ax1.scatter(iq_i, iq_q, c=colors, s=40, edgecolors="k",
                    linewidths=0.4, alpha=0.7)
        ax1.axhline(0, color="k", linewidth=0.5)
        ax1.axvline(0, color="k", linewidth=0.5)
        ax1.set_title("Constellation Diagram")
        ax1.set_xlabel("In-phase (I)")
        ax1.set_ylabel("Quadrature (Q)")
        ax1.set_aspect("equal")
        ax1.grid(True, alpha=0.3)

        # Time-domain snippet (first 8 bits)
        show = min(8, num_bits) * sps
        t_ms = np.arange(show) / fs * 1000
        ax2.plot(t_ms, noisy[:show], color="tab:blue", linewidth=0.6,
                 label="Received (noisy)")
        ax2.plot(t_ms, modulated[:show], color="tab:orange", linewidth=1.0,
                 alpha=0.6, label="Clean signal")
        ax2.set_title("Waveform (first 8 bits)")
        ax2.set_xlabel("Time (ms)")
        ax2.set_ylabel("Amplitude")
        ax2.legend(fontsize=8)
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

    interact(
        _update,
        snr_db=FloatSlider(value=20, min=-5, max=40, step=1,
                           description="SNR (dB)"),
        num_bits=IntSlider(value=50, min=10, max=500, step=10,
                           description="Num bits"),
    )


# ---------------------------------------------------------------------------
# Standalone execution
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    matplotlib.use("Agg")
    fig = plot_bpsk_demo()
    fig.savefig("bpsk_demo.png", dpi=150)
    print("Saved bpsk_demo.png")
