"""
Quadrature Amplitude Modulation (QAM)

QAM encodes digital data by varying BOTH the amplitude AND phase of a carrier
signal. Think of it as combining two ideas:
  - ASK (Amplitude Shift Keying): encode data by changing signal strength
  - PSK (Phase Shift Keying): encode data by changing signal timing/angle

Each transmitted "symbol" is a unique combination of amplitude and phase,
represented as a point on a 2D grid called a "constellation":
  - 4-QAM  (= QPSK):  4 points,  2 bits per symbol
  - 16-QAM:           16 points,  4 bits per symbol (4x4 grid)
  - 64-QAM:           64 points,  6 bits per symbol (8x8 grid)
  - 256-QAM:         256 points,  8 bits per symbol (16x16 grid)

The two axes of the constellation are called I (In-phase) and Q (Quadrature).
These correspond to two carrier signals 90 degrees apart — one is a cosine,
the other a sine. Each symbol sets the amplitude of both carriers, and the
receiver recovers the data by measuring where the received signal lands on
the I/Q grid.

Higher-order QAM packs more bits into each symbol (faster data rate), but the
constellation points are closer together, making the signal more vulnerable
to noise. This is the fundamental tradeoff in digital communications.

Gray coding arranges the bit labels so that adjacent constellation points
differ by only 1 bit. This way, when noise pushes a received symbol to a
neighboring point, only 1 bit is wrong instead of potentially several.

Real-world uses: WiFi (802.11ac/ax uses up to 1024-QAM), 4G/5G cellular,
cable TV (QAM-256), and digital satellite.
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Core functions (importable by notebook)
# ---------------------------------------------------------------------------

def _gray_code(n):
    """Return Gray code for integer n: adjacent values differ by 1 bit."""
    return n ^ (n >> 1)


def _gray_code_table(bits_per_dim):
    """Build a Gray-coded ordering for one dimension of the constellation.

    For a 16-QAM constellation (4x4 grid), bits_per_dim = 2, so the axis
    indices [0, 1, 2, 3] are reordered by Gray code: [0, 1, 3, 2].
    """
    size = 2 ** bits_per_dim
    return [_gray_code(i) for i in range(size)]


def qam_constellation(M):
    """Generate the ideal QAM constellation with Gray-coded bit labels.

    Parameters
    ----------
    M : int
        Constellation order — must be a power of 4 (4, 16, 64, 256, ...).

    Returns
    -------
    symbols : ndarray of complex, shape (M,)
        Constellation points (I + jQ), normalized to unit average power.
    bit_labels : list of str
        Gray-coded bit string for each constellation point.
    """
    k = int(np.log2(M))
    if 2 ** k != M or k % 2 != 0:
        raise ValueError(f"M must be a power of 4 (4, 16, 64, ...), got {M}")

    bits_per_dim = k // 2
    n_side = int(np.sqrt(M))

    # PAM levels centered at zero: e.g. for 4x4 -> [-3, -1, 1, 3]
    levels = np.arange(n_side) * 2 - (n_side - 1)

    # Gray-coded index ordering for each axis
    gray_order = _gray_code_table(bits_per_dim)

    symbols = np.empty(M, dtype=complex)
    bit_labels = [''] * M

    idx = 0
    for i_gray_idx, i_level_idx in enumerate(gray_order):
        for q_gray_idx, q_level_idx in enumerate(gray_order):
            i_val = levels[i_level_idx]
            q_val = levels[q_level_idx]
            symbols[idx] = i_val + 1j * q_val

            # Bit label: MSBs select I level, LSBs select Q level
            i_bits = format(i_gray_idx, f'0{bits_per_dim}b')
            q_bits = format(q_gray_idx, f'0{bits_per_dim}b')
            bit_labels[idx] = i_bits + q_bits
            idx += 1

    # Normalize to unit average power
    avg_power = np.mean(np.abs(symbols) ** 2)
    symbols /= np.sqrt(avg_power)

    return symbols, bit_labels


def qam_modulate(bits, M, fs=1.0, sps=8):
    """Modulate a bit stream using M-QAM.

    Parameters
    ----------
    bits : array-like of int
        Binary data (0s and 1s). Length must be a multiple of log2(M).
    M : int
        QAM order (4, 16, 64, 256).
    fs : float
        Sample rate in samples per second (default 1.0 for baseband).
    sps : int
        Samples per symbol — how many output samples represent each symbol.
        Higher values give a smoother waveform.

    Returns
    -------
    signal : ndarray of complex
        Modulated baseband signal (I + jQ samples).
    symbols : ndarray of complex
        The constellation points chosen for each symbol period.
    """
    bits = np.array(bits, dtype=int)
    k = int(np.log2(M))

    if len(bits) % k != 0:
        raise ValueError(
            f"Number of bits ({len(bits)}) must be a multiple of "
            f"log2(M) = {k}"
        )

    constellation, labels = qam_constellation(M)

    # Build a lookup: bit string -> constellation point
    label_to_symbol = {lbl: sym for lbl, sym in zip(labels, constellation)}

    # Group bits into k-bit chunks, map each to a constellation point
    num_symbols = len(bits) // k
    symbols = np.empty(num_symbols, dtype=complex)
    for i in range(num_symbols):
        chunk = bits[i * k:(i + 1) * k]
        bit_str = ''.join(str(b) for b in chunk)
        symbols[i] = label_to_symbol[bit_str]

    # Upsample: insert each symbol followed by (sps-1) zeros, then
    # apply a simple rectangular pulse (repeat each symbol sps times).
    signal = np.repeat(symbols, sps)

    return signal, symbols


def qam_demodulate(signal, M, fs=1.0, sps=8):
    """Demodulate an M-QAM baseband signal back to bits.

    Uses minimum-distance decoding: each received sample (taken once per
    symbol period) is matched to the nearest ideal constellation point.

    Parameters
    ----------
    signal : ndarray of complex
        Received baseband signal.
    M : int
        QAM order.
    fs : float
        Sample rate (must match modulator).
    sps : int
        Samples per symbol (must match modulator).

    Returns
    -------
    bits : ndarray of int
        Recovered binary data.
    received_symbols : ndarray of complex
        The sample values at each symbol instant (for plotting).
    """
    constellation, labels = qam_constellation(M)
    k = int(np.log2(M))

    # Sample at the center of each symbol period
    sample_indices = np.arange(sps // 2, len(signal), sps)
    received_symbols = signal[sample_indices]

    # Minimum-distance decoding
    bits = []
    for sym in received_symbols:
        distances = np.abs(constellation - sym)
        nearest_idx = np.argmin(distances)
        bits.extend(int(b) for b in labels[nearest_idx])

    return np.array(bits, dtype=int), received_symbols


def add_awgn(signal, snr_db):
    """Add Additive White Gaussian Noise to a complex signal.

    AWGN is the simplest noise model: random noise added uniformly across
    all frequencies. "White" means equal power at all frequencies (like
    white light contains all colors). "Gaussian" means the noise amplitude
    follows a bell curve.

    Parameters
    ----------
    signal : ndarray of complex
        Clean input signal.
    snr_db : float
        Signal-to-noise ratio in decibels. Higher = less noise.
        Typical values: 30 dB (very clean), 20 dB (good), 10 dB (noisy).

    Returns
    -------
    noisy_signal : ndarray of complex
        Signal with noise added.
    """
    signal_power = np.mean(np.abs(signal) ** 2)
    snr_linear = 10 ** (snr_db / 10)
    noise_power = signal_power / snr_linear
    noise = np.sqrt(noise_power / 2) * (
        np.random.randn(len(signal)) + 1j * np.random.randn(len(signal))
    )
    return signal + noise


def compute_ber(M, snr_db, num_bits=100000):
    """Compute Bit Error Rate for M-QAM at a given SNR.

    Transmits random bits through an AWGN channel and counts errors.

    Parameters
    ----------
    M : int
        QAM order.
    snr_db : float
        Channel SNR in dB.
    num_bits : int
        Number of bits to simulate (more = smoother BER curve).

    Returns
    -------
    ber : float
        Measured bit error rate (0.0 = perfect, 0.5 = useless).
    """
    k = int(np.log2(M))
    # Round down to multiple of k
    num_bits = (num_bits // k) * k

    bits_tx = np.random.randint(0, 2, num_bits)
    signal, _ = qam_modulate(bits_tx, M, sps=1)
    noisy = add_awgn(signal, snr_db)
    bits_rx, _ = qam_demodulate(noisy, M, sps=1)

    min_len = min(len(bits_tx), len(bits_rx))
    errors = np.sum(bits_tx[:min_len] != bits_rx[:min_len])
    return errors / min_len


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def plot_qam_demo():
    """Generate a multi-panel figure showing QAM concepts.

    Panel layout:
      Top-left:     16-QAM constellation (clean) with Gray-coded bit labels
      Top-right:    16-QAM constellation at various SNR levels
      Bottom-left:  64-QAM constellation (clean)
      Bottom-right: BER vs SNR comparison for 4, 16, and 64-QAM
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle("Quadrature Amplitude Modulation (QAM)", fontsize=15)

    # --- Panel 1: Clean 16-QAM with bit labels ---
    ax = axes[0, 0]
    symbols_16, labels_16 = qam_constellation(16)
    ax.scatter(symbols_16.real, symbols_16.imag, s=100, c="tab:blue",
               zorder=5, edgecolors="black", linewidths=0.5)
    for sym, lbl in zip(symbols_16, labels_16):
        ax.annotate(lbl, (sym.real, sym.imag), textcoords="offset points",
                    xytext=(0, 10), ha="center", fontsize=7,
                    fontfamily="monospace")
    ax.set_title("16-QAM Constellation (Gray Coded)")
    ax.set_xlabel("In-phase (I)")
    ax.set_ylabel("Quadrature (Q)")
    ax.axhline(0, color="gray", linewidth=0.5)
    ax.axvline(0, color="gray", linewidth=0.5)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)

    # --- Panel 2: 16-QAM with noise at various SNR levels ---
    ax = axes[0, 1]
    snr_levels = [25, 15, 8]
    colors = ["tab:green", "tab:orange", "tab:red"]
    num_syms = 500
    bits_demo = np.random.randint(0, 2, num_syms * 4)

    for snr, color in zip(snr_levels, colors):
        sig, _ = qam_modulate(bits_demo, 16, sps=1)
        noisy = add_awgn(sig, snr)
        ax.scatter(noisy.real, noisy.imag, s=8, alpha=0.4, c=color,
                   label=f"SNR = {snr} dB")

    # Overlay ideal points
    ax.scatter(symbols_16.real, symbols_16.imag, s=80, c="black",
               marker="x", zorder=10, linewidths=2, label="Ideal")
    ax.set_title("16-QAM with AWGN Noise")
    ax.set_xlabel("In-phase (I)")
    ax.set_ylabel("Quadrature (Q)")
    ax.axhline(0, color="gray", linewidth=0.5)
    ax.axvline(0, color="gray", linewidth=0.5)
    ax.set_aspect("equal")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, alpha=0.3)

    # --- Panel 3: 64-QAM constellation ---
    ax = axes[1, 0]
    symbols_64, labels_64 = qam_constellation(64)
    ax.scatter(symbols_64.real, symbols_64.imag, s=40, c="tab:purple",
               zorder=5, edgecolors="black", linewidths=0.3)
    for sym, lbl in zip(symbols_64, labels_64):
        ax.annotate(lbl, (sym.real, sym.imag), textcoords="offset points",
                    xytext=(0, 7), ha="center", fontsize=4.5,
                    fontfamily="monospace", color="dimgray")
    ax.set_title("64-QAM Constellation (8x8 grid, 6 bits/symbol)")
    ax.set_xlabel("In-phase (I)")
    ax.set_ylabel("Quadrature (Q)")
    ax.axhline(0, color="gray", linewidth=0.5)
    ax.axvline(0, color="gray", linewidth=0.5)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)

    # --- Panel 4: BER vs SNR ---
    ax = axes[1, 1]
    snr_range = np.arange(0, 28, 1.5)
    qam_orders = [4, 16, 64]
    markers = ["o", "s", "^"]
    colors_ber = ["tab:blue", "tab:orange", "tab:green"]

    for M, marker, color in zip(qam_orders, markers, colors_ber):
        ber_values = []
        for snr in snr_range:
            ber = compute_ber(M, snr, num_bits=50000)
            ber_values.append(max(ber, 1e-6))  # floor for log plot
        ax.semilogy(snr_range, ber_values, marker=marker, markersize=4,
                    color=color, label=f"{M}-QAM ({int(np.log2(M))} bits/sym)")

    ax.set_title("Bit Error Rate vs SNR")
    ax.set_xlabel("SNR (dB)")
    ax.set_ylabel("Bit Error Rate (BER)")
    ax.set_ylim(1e-6, 1)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, which="both")
    ax.axhline(1e-3, color="red", linestyle="--", alpha=0.4, linewidth=0.8)
    ax.annotate("Typical target BER", xy=(1, 1e-3), fontsize=7,
                color="red", alpha=0.6)

    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Interactive widget (for Jupyter notebook)
# ---------------------------------------------------------------------------

def interactive_qam():
    """Launch ipywidgets sliders for exploring QAM constellations.

    Sliders:
      - QAM Order: choose between 4, 16, or 64-QAM
      - SNR (dB): control noise level from 0 (very noisy) to 30 (clean)

    The display shows the noisy constellation with ideal points overlaid,
    and reports the measured bit error rate.
    """
    try:
        from ipywidgets import interact, IntSlider, FloatSlider
    except ImportError:
        print("ipywidgets not available — run in Jupyter for interactive mode")
        return

    def _update(order_index=1, snr_db=20.0):
        order_map = {0: 4, 1: 16, 2: 64}
        M = order_map[order_index]
        k = int(np.log2(M))
        num_symbols = 600
        num_bits = num_symbols * k

        bits_tx = np.random.randint(0, 2, num_bits)
        signal, _ = qam_modulate(bits_tx, M, sps=1)
        noisy = add_awgn(signal, snr_db)
        bits_rx, rx_symbols = qam_demodulate(noisy, M, sps=1)

        min_len = min(len(bits_tx), len(bits_rx))
        ber = np.sum(bits_tx[:min_len] != bits_rx[:min_len]) / min_len

        constellation, labels = qam_constellation(M)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

        # Left: noisy constellation
        ax1.scatter(rx_symbols.real, rx_symbols.imag, s=12, alpha=0.5,
                    c="tab:blue", label="Received")
        ax1.scatter(constellation.real, constellation.imag, s=80, c="red",
                    marker="x", zorder=10, linewidths=2, label="Ideal")
        ax1.set_title(f"{M}-QAM  |  SNR = {snr_db:.0f} dB  |  "
                      f"BER = {ber:.4f}")
        ax1.set_xlabel("In-phase (I)")
        ax1.set_ylabel("Quadrature (Q)")
        ax1.axhline(0, color="gray", linewidth=0.5)
        ax1.axvline(0, color="gray", linewidth=0.5)
        ax1.set_aspect("equal")
        ax1.legend(fontsize=8)
        ax1.grid(True, alpha=0.3)

        # Right: ideal constellation with bit labels
        ax2.scatter(constellation.real, constellation.imag, s=100,
                    c="tab:green", zorder=5, edgecolors="black",
                    linewidths=0.5)
        for sym, lbl in zip(constellation, labels):
            ax2.annotate(lbl, (sym.real, sym.imag),
                         textcoords="offset points", xytext=(0, 10),
                         ha="center", fontsize=7 if M <= 16 else 5,
                         fontfamily="monospace")
        ax2.set_title(f"{M}-QAM Ideal Constellation "
                      f"({k} bits/symbol, {M} points)")
        ax2.set_xlabel("In-phase (I)")
        ax2.set_ylabel("Quadrature (Q)")
        ax2.axhline(0, color="gray", linewidth=0.5)
        ax2.axvline(0, color="gray", linewidth=0.5)
        ax2.set_aspect("equal")
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

    interact(
        _update,
        order_index=IntSlider(value=1, min=0, max=2, step=1,
                              description="QAM Order",
                              readout=False,
                              style={"description_width": "initial"}),
        snr_db=FloatSlider(value=20.0, min=0.0, max=30.0, step=1.0,
                           description="SNR (dB)"),
    )


# ---------------------------------------------------------------------------
# Standalone execution
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    matplotlib.use("Agg")

    fig = plot_qam_demo()
    fig.savefig("qam_demo.png", dpi=150)
    print("Saved qam_demo.png")

    plt.show()
