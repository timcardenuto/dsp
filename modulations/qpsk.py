"""
QPSK — Quadrature Phase Shift Keying

BPSK encodes one bit per symbol by flipping a carrier's phase between 0 and
180 degrees. QPSK doubles the data rate by using FOUR phase states (45, 135,
225, 315 degrees), encoding TWO bits per symbol. The trick is splitting the
carrier into two independent channels:

  - I (In-phase):     the cosine component — carries one bit
  - Q (Quadrature):   the sine component  — carries the other bit

Each channel is essentially its own BPSK signal. Because cosine and sine are
orthogonal (their product integrates to zero over a full cycle), the two
channels don't interfere with each other. The receiver can separate them
cleanly.

The result: QPSK achieves TWICE the spectral efficiency of BPSK — same
bandwidth, double the data rate — at the cost of a slightly more complex
transmitter and receiver.

Constellation diagram
---------------------
The four symbols sit at equal distances from the origin, forming a square:

        Q
        |
   01   |   11        Each point = one 2-bit symbol
   -----+----->  I    Distance from origin = signal amplitude
   00   |   10        Angle = carrier phase
        |

Gray coding (adjacent symbols differ by only 1 bit) minimizes bit errors
when noise pushes a received symbol to a neighboring point.

Converted from QPSK.m (MATLAB reference implementation).
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Core functions (importable by notebook)
# ---------------------------------------------------------------------------

# Gray-coded symbol map: index = symbol (0-3), value = (I, Q)
# Bit pair -> symbol:  00->0, 01->1, 11->2, 10->3
# Phase angles: 225, 135, 45, 315 degrees
_GRAY_MAP = {
    (0, 0): (-1, -1),   # 225 degrees
    (0, 1): (-1, +1),   # 135 degrees
    (1, 1): (+1, +1),   #  45 degrees
    (1, 0): (+1, -1),   # 315 degrees
}

# Reverse map for demodulation: (I sign, Q sign) -> (bit0, bit1)
_GRAY_DEMAP = {v: k for k, v in _GRAY_MAP.items()}


def qpsk_modulate(bits, carrier_freq, fs, sps):
    """Modulate a bit sequence onto a QPSK carrier.

    Takes raw bits, groups them into pairs, and produces a real-valued
    passband signal: each symbol interval contains a carrier whose phase
    encodes two bits.

    Parameters
    ----------
    bits : array-like – binary data (0s and 1s). Length must be even.
    carrier_freq : float – carrier frequency in Hz
    fs : float – sample rate in samples/second
    sps : int – samples per symbol (controls symbol duration)

    Returns
    -------
    t : ndarray – time vector for the full signal
    signal : ndarray – modulated QPSK waveform
    i_component : ndarray – in-phase (cosine) component alone
    q_component : ndarray – quadrature (sine) component alone
    symbols : list of tuples – (I, Q) values per symbol for constellation
    """
    bits = np.array(bits, dtype=int)
    if len(bits) % 2 != 0:
        raise ValueError("Bit sequence length must be even for QPSK")

    num_symbols = len(bits) // 2
    t_symbol = np.arange(sps) / fs            # time vector for one symbol
    t = np.arange(num_symbols * sps) / fs      # full time vector

    signal = np.zeros(len(t))
    i_component = np.zeros(len(t))
    q_component = np.zeros(len(t))
    symbols = []

    for k in range(num_symbols):
        b0, b1 = bits[2 * k], bits[2 * k + 1]
        i_val, q_val = _GRAY_MAP[(b0, b1)]
        symbols.append((i_val, q_val))

        idx_start = k * sps
        idx_end = idx_start + sps

        # In-phase: data * cos(2*pi*fc*t)
        i_wave = i_val * np.cos(2 * np.pi * carrier_freq * t_symbol)
        # Quadrature: data * sin(2*pi*fc*t)
        q_wave = q_val * np.sin(2 * np.pi * carrier_freq * t_symbol)

        i_component[idx_start:idx_end] = i_wave
        q_component[idx_start:idx_end] = q_wave
        signal[idx_start:idx_end] = i_wave + q_wave

    return t, signal, i_component, q_component, symbols


def qpsk_demodulate(signal, carrier_freq, fs, sps):
    """Demodulate a QPSK signal back to bits.

    Multiplies the received signal by cosine and sine references, integrates
    over each symbol period (matched filter), and decides which 2-bit symbol
    was sent based on the sign of each integral.

    Parameters
    ----------
    signal : ndarray – received QPSK waveform
    carrier_freq : float – carrier frequency in Hz (must match transmitter)
    fs : float – sample rate in samples/second
    sps : int – samples per symbol

    Returns
    -------
    bits : ndarray – recovered bit sequence
    symbols : list of tuples – detected (I, Q) values per symbol
    """
    num_symbols = len(signal) // sps
    t_symbol = np.arange(sps) / fs

    bits = []
    symbols = []

    for k in range(num_symbols):
        idx_start = k * sps
        idx_end = idx_start + sps
        chunk = signal[idx_start:idx_end]

        # Correlate with cosine and sine references
        i_corr = np.sum(chunk * np.cos(2 * np.pi * carrier_freq * t_symbol))
        q_corr = np.sum(chunk * np.sin(2 * np.pi * carrier_freq * t_symbol))

        # Hard decision: sign determines the bit
        i_decision = +1 if i_corr >= 0 else -1
        q_decision = +1 if q_corr >= 0 else -1

        symbols.append((i_decision, q_decision))
        b0, b1 = _GRAY_DEMAP[(i_decision, q_decision)]
        bits.extend([b0, b1])

    return np.array(bits, dtype=int), symbols


def add_awgn(signal, snr_db):
    """Add white Gaussian noise to achieve a target SNR (in dB).

    SNR is the ratio of signal power to noise power. At 0 dB the noise is
    as strong as the signal; at 20 dB the signal is 100x stronger than noise.

    Parameters
    ----------
    signal : ndarray – clean signal
    snr_db : float – desired signal-to-noise ratio in decibels

    Returns
    -------
    noisy : ndarray – signal with added noise
    noise : ndarray – the noise that was added (useful for analysis)
    """
    sig_power = np.mean(signal ** 2)
    noise_power = sig_power / (10 ** (snr_db / 10))
    noise = np.sqrt(noise_power) * np.random.randn(len(signal))
    return signal + noise, noise


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def plot_qpsk_demo(carrier_freq=5e3, fs=200e3, sps=200, snr_db=10):
    """Generate a multi-panel QPSK overview figure.

    Panels
    ------
    1. I (in-phase) component waveform
    2. Q (quadrature) component waveform
    3. Combined QPSK signal
    4. Constellation diagram — clean vs noisy
    5. Frequency spectrum of the modulated signal
    """
    # Example data: 10 bits -> 5 QPSK symbols
    bits = np.array([0, 1, 0, 1, 1, 1, 0, 0, 1, 1])

    t, signal, i_comp, q_comp, symbols = qpsk_modulate(
        bits, carrier_freq, fs, sps
    )
    t_ms = t * 1000  # convert to milliseconds for display

    # Add noise for constellation comparison
    noisy_signal, _ = add_awgn(signal, snr_db)
    _, noisy_symbols = qpsk_demodulate(noisy_signal, carrier_freq, fs, sps)

    # Compute spectrum
    N = len(signal)
    Y = np.fft.fft(signal)
    freqs = np.fft.fftfreq(N, 1 / fs)[:N // 2]
    magnitude = 2.0 / N * np.abs(Y[:N // 2])

    # --- Build figure ---
    fig, axes = plt.subplots(3, 2, figsize=(14, 10))
    fig.suptitle("QPSK Modulation Overview", fontsize=14)

    # Panel 1: I component
    ax = axes[0, 0]
    ax.plot(t_ms, i_comp, color="tab:blue")
    ax.set_title("I (In-Phase) Component")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude")
    ax.grid(True, alpha=0.3)
    # Annotate bit pairs
    for k in range(len(bits) // 2):
        mid = (k + 0.5) * sps / fs * 1000
        ax.text(mid, 1.15, f"{bits[2*k]}{bits[2*k+1]}",
                ha="center", fontsize=8, color="tab:blue")

    # Panel 2: Q component
    ax = axes[0, 1]
    ax.plot(t_ms, q_comp, color="tab:orange")
    ax.set_title("Q (Quadrature) Component")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude")
    ax.grid(True, alpha=0.3)

    # Panel 3: Combined QPSK signal
    ax = axes[1, 0]
    ax.plot(t_ms, signal, color="tab:green")
    ax.set_title("Combined QPSK Signal (I + Q)")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude")
    ax.grid(True, alpha=0.3)

    # Panel 4: Combined QPSK signal with noise
    ax = axes[1, 1]
    ax.plot(t_ms, noisy_signal, color="tab:red", alpha=0.7)
    ax.set_title(f"QPSK Signal + Noise (SNR = {snr_db} dB)")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude")
    ax.grid(True, alpha=0.3)

    # Panel 5: Constellation diagram
    ax = axes[2, 0]
    # Clean constellation
    sym_i = [s[0] for s in symbols]
    sym_q = [s[1] for s in symbols]
    ax.scatter(sym_i, sym_q, s=120, c="tab:blue", marker="o",
               zorder=3, label="Clean")
    # Noisy constellation — generate many symbols for visual spread
    np.random.seed(42)
    many_bits = np.random.randint(0, 2, 200)
    _, many_sig, _, _, _ = qpsk_modulate(many_bits, carrier_freq, fs, sps)
    noisy_many, _ = add_awgn(many_sig, snr_db)
    _, noisy_syms = qpsk_demodulate(noisy_many, carrier_freq, fs, sps)
    # Recover continuous I/Q correlation values (not hard decisions) for scatter
    num_syms = len(many_sig) // sps
    t_sym = np.arange(sps) / fs
    soft_i, soft_q = [], []
    for k in range(num_syms):
        chunk = noisy_many[k * sps:(k + 1) * sps]
        soft_i.append(np.sum(chunk * np.cos(2 * np.pi * carrier_freq * t_sym)))
        soft_q.append(np.sum(chunk * np.sin(2 * np.pi * carrier_freq * t_sym)))
    # Normalize for display
    scale = sps / 2
    soft_i = np.array(soft_i) / scale
    soft_q = np.array(soft_q) / scale
    ax.scatter(soft_i, soft_q, s=15, c="tab:red", alpha=0.4, label="Noisy")
    ax.axhline(0, color="gray", linewidth=0.5)
    ax.axvline(0, color="gray", linewidth=0.5)
    ax.set_title("Constellation Diagram")
    ax.set_xlabel("In-Phase (I)")
    ax.set_ylabel("Quadrature (Q)")
    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim(-2.5, 2.5)
    ax.set_aspect("equal")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    # Label ideal points
    for (b0, b1), (iv, qv) in _GRAY_MAP.items():
        ax.annotate(f"{b0}{b1}", (iv, qv), textcoords="offset points",
                    xytext=(10, 10), fontsize=9, fontweight="bold")

    # Panel 6: Frequency spectrum
    ax = axes[2, 1]
    ax.plot(freqs / 1000, magnitude, color="tab:purple")
    ax.set_title("Frequency Spectrum")
    ax.set_xlabel("Frequency (kHz)")
    ax.set_ylabel("Magnitude")
    ax.set_xlim(0, carrier_freq * 3 / 1000)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Interactive widget (for Jupyter notebook)
# ---------------------------------------------------------------------------

def interactive_qpsk():
    """Launch ipywidgets sliders for exploring QPSK under varying noise.

    Sliders
    -------
    SNR (dB) : controls how much noise is added to the signal. Lower values
               make the constellation points spread out and eventually cause
               bit errors.
    """
    try:
        from ipywidgets import interact, FloatSlider
    except ImportError:
        print("ipywidgets not available — run in Jupyter for interactive mode")
        return

    def _update(snr_db=15.0):
        carrier_freq = 5e3
        fs = 200e3
        sps = 200

        np.random.seed(0)
        bits = np.random.randint(0, 2, 200)
        t, signal, i_comp, q_comp, symbols = qpsk_modulate(
            bits, carrier_freq, fs, sps
        )

        noisy_signal, _ = add_awgn(signal, snr_db)
        recovered_bits, _ = qpsk_demodulate(noisy_signal, carrier_freq, fs, sps)

        ber = np.mean(bits != recovered_bits)

        # Soft I/Q values for constellation scatter
        num_syms = len(signal) // sps
        t_sym = np.arange(sps) / fs
        soft_i, soft_q = [], []
        for k in range(num_syms):
            chunk = noisy_signal[k * sps:(k + 1) * sps]
            soft_i.append(
                np.sum(chunk * np.cos(2 * np.pi * carrier_freq * t_sym))
            )
            soft_q.append(
                np.sum(chunk * np.sin(2 * np.pi * carrier_freq * t_sym))
            )
        scale = sps / 2
        soft_i = np.array(soft_i) / scale
        soft_q = np.array(soft_q) / scale

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        # Constellation
        ax1.scatter(soft_i, soft_q, s=20, c="tab:red", alpha=0.5)
        # Ideal points
        for (b0, b1), (iv, qv) in _GRAY_MAP.items():
            ax1.scatter(iv, qv, s=100, c="tab:blue", zorder=3)
            ax1.annotate(f"{b0}{b1}", (iv, qv), textcoords="offset points",
                         xytext=(10, 10), fontsize=9, fontweight="bold")
        ax1.axhline(0, color="gray", linewidth=0.5)
        ax1.axvline(0, color="gray", linewidth=0.5)
        ax1.set_title(f"Constellation (SNR = {snr_db:.0f} dB, BER = {ber:.3f})")
        ax1.set_xlabel("In-Phase (I)")
        ax1.set_ylabel("Quadrature (Q)")
        ax1.set_xlim(-3, 3)
        ax1.set_ylim(-3, 3)
        ax1.set_aspect("equal")
        ax1.grid(True, alpha=0.3)

        # Time domain comparison
        show_samples = min(2000, len(signal))
        t_ms = t[:show_samples] * 1000
        ax2.plot(t_ms, signal[:show_samples], alpha=0.7, label="Clean")
        ax2.plot(t_ms, noisy_signal[:show_samples], alpha=0.5, label="Noisy")
        ax2.set_title("QPSK Waveform")
        ax2.set_xlabel("Time (ms)")
        ax2.set_ylabel("Amplitude")
        ax2.legend(fontsize=8)
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

    interact(
        _update,
        snr_db=FloatSlider(
            value=15.0, min=-5.0, max=30.0, step=1.0,
            description="SNR (dB)",
            style={"description_width": "initial"},
        ),
    )


# ---------------------------------------------------------------------------
# Standalone execution
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    matplotlib.use("Agg")

    fig = plot_qpsk_demo()
    fig.savefig("qpsk_modulation.png", dpi=150)
    print("Saved qpsk_modulation.png")

    plt.show()
