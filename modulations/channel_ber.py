"""
Channel Effects & Bit Error Rate (BER)

When you send digital data wirelessly (or over any medium), the signal travels
through a "channel" — air, copper cable, fiber optic, etc.  The channel is
never perfect: it adds noise, distortion, and interference.  The central
question of digital communications is: given these impairments, how many bits
arrive correctly?

Key concepts for non-RF people
-------------------------------
AWGN (Additive White Gaussian Noise):
    The simplest and most common noise model.  "Additive" means the noise is
    added on top of your signal.  "White" means it has equal power at all
    frequencies (like white light contains all colors).  "Gaussian" means the
    amplitude of each noise sample follows a bell curve.  Every real system has
    at least this much noise — it comes from thermal motion of electrons.

SNR (Signal-to-Noise Ratio):
    How much stronger your signal is compared to the noise, measured in
    decibels (dB).  At 0 dB the signal and noise are equally strong; at 20 dB
    the signal is 100x more powerful.  Higher SNR = cleaner reception.

BER (Bit Error Rate):
    The fraction of bits that arrive incorrectly.  A BER of 1e-3 means one bit
    in a thousand is wrong.  This is THE single most important metric for any
    digital communication link.  As SNR decreases, BER increases (more errors).

Modulation trade-offs:
    Different modulation schemes pack different numbers of bits into each
    symbol.  BPSK sends 1 bit/symbol and is very robust; 64-QAM sends 6
    bits/symbol but needs a much cleaner channel.  The BER-vs-SNR curve shows
    this trade-off clearly: higher-order modulations need more SNR to achieve
    the same error rate.

Constellation degradation:
    A constellation diagram shows where received symbols land on the I/Q plane.
    With no noise, they sit on perfect grid points.  As noise increases, the
    points spread into clouds.  When a cloud overlaps the decision boundary,
    the receiver picks the wrong symbol — that is a bit error.

Converted from WiFi_sim.m (MATLAB reference implementation).
"""

import numpy as np
from scipy.special import erfc
import matplotlib
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Core functions (importable by notebook)
# ---------------------------------------------------------------------------

def add_awgn(signal, snr_db):
    """Add white Gaussian noise to achieve a target SNR (in dB).

    Think of it as layering static on top of your signal.  The noise power is
    calculated so the ratio of signal power to noise power equals the requested
    SNR.

    Parameters
    ----------
    signal : ndarray – clean (complex or real) signal samples
    snr_db : float – desired signal-to-noise ratio in decibels

    Returns
    -------
    noisy : ndarray – signal with noise added
    noise : ndarray – the noise that was added (useful for analysis)
    """
    signal = np.asarray(signal)
    sig_power = np.mean(np.abs(signal) ** 2)
    noise_power = sig_power / (10 ** (snr_db / 10))
    # For complex signals, split noise equally between real and imaginary
    if np.iscomplexobj(signal):
        noise = np.sqrt(noise_power / 2) * (
            np.random.randn(len(signal)) + 1j * np.random.randn(len(signal))
        )
    else:
        noise = np.sqrt(noise_power) * np.random.randn(len(signal))
    return signal + noise, noise


def compute_ber(tx_bits, rx_bits):
    """Compute the Bit Error Rate between transmitted and received bits.

    Simply counts how many bits differ and divides by the total number of bits.
    A BER of 0 means perfect reception; 0.5 means you might as well flip coins.

    Parameters
    ----------
    tx_bits : array-like – original transmitted bits (0s and 1s)
    rx_bits : array-like – received/decoded bits (0s and 1s)

    Returns
    -------
    ber : float – fraction of bits in error (0.0 to 1.0)
    num_errors : int – total count of bit errors
    """
    tx = np.asarray(tx_bits, dtype=int)
    rx = np.asarray(rx_bits, dtype=int)
    if len(tx) != len(rx):
        raise ValueError(
            f"Bit sequences must be same length (got {len(tx)} vs {len(rx)})"
        )
    num_errors = int(np.sum(tx != rx))
    ber = num_errors / len(tx)
    return ber, num_errors


def bpsk_theoretical_ber(snr_db):
    """Theoretical BER for BPSK over an AWGN channel.

    This is the closed-form "best you can do" formula from communication
    theory.  It uses the complementary error function (erfc).

    Parameters
    ----------
    snr_db : float or ndarray – SNR values in dB

    Returns
    -------
    ber : float or ndarray – theoretical bit error rate(s)
    """
    snr_linear = 10 ** (np.asarray(snr_db, dtype=float) / 10)
    return 0.5 * erfc(np.sqrt(snr_linear))


def qpsk_theoretical_ber(snr_db):
    """Theoretical BER for QPSK over an AWGN channel.

    QPSK has the same BER-per-bit as BPSK (because the I and Q channels are
    independent BPSK signals), but carries 2 bits per symbol.  The formula
    is identical to BPSK when expressed in terms of Eb/No (energy per bit).

    Parameters
    ----------
    snr_db : float or ndarray – SNR per bit in dB

    Returns
    -------
    ber : float or ndarray – theoretical bit error rate(s)
    """
    # For QPSK with Gray coding, BER = BPSK BER
    return bpsk_theoretical_ber(snr_db)


def qam_theoretical_ber(snr_db, M):
    """Approximate theoretical BER for M-QAM over an AWGN channel.

    Uses the standard approximation valid for square QAM constellations
    (M = 4, 16, 64, 256, ...) with Gray coding.

    Parameters
    ----------
    snr_db : float or ndarray – SNR per bit in dB
    M : int – constellation size (must be a power of 4)

    Returns
    -------
    ber : float or ndarray – approximate theoretical bit error rate(s)
    """
    k = np.log2(M)  # bits per symbol
    snr_linear = 10 ** (np.asarray(snr_db, dtype=float) / 10)
    # Convert Eb/No to Es/No (energy per symbol to noise)
    es_no = snr_linear * k
    ber = (2 * (1 - 1 / np.sqrt(M)) / k) * erfc(
        np.sqrt(3 * es_no / (2 * (M - 1)))
    )
    return np.clip(ber, 0, 0.5)


# ---------------------------------------------------------------------------
# Baseband modulation / demodulation helpers (complex I/Q, no carrier)
# ---------------------------------------------------------------------------

def _bpsk_mod(bits):
    """BPSK: 0 -> -1, 1 -> +1 (real-valued)."""
    return 2.0 * np.asarray(bits, dtype=float) - 1.0


def _bpsk_demod(symbols):
    """BPSK hard decision: positive -> 1, negative -> 0."""
    return (np.real(symbols) >= 0).astype(int)


def _qpsk_mod(bits):
    """QPSK: pairs of bits -> complex symbols on unit circle."""
    bits = np.asarray(bits, dtype=int)
    if len(bits) % 2 != 0:
        bits = np.append(bits, 0)  # pad
    I = 2.0 * bits[0::2] - 1.0
    Q = 2.0 * bits[1::2] - 1.0
    return (I + 1j * Q) / np.sqrt(2)


def _qpsk_demod(symbols):
    """QPSK hard decision: sign of I and Q components."""
    I_bits = (np.real(symbols) >= 0).astype(int)
    Q_bits = (np.imag(symbols) >= 0).astype(int)
    bits = np.empty(2 * len(symbols), dtype=int)
    bits[0::2] = I_bits
    bits[1::2] = Q_bits
    return bits


def _qam_constellation(M):
    """Generate a normalized M-QAM constellation (Gray-coded where practical).

    Returns an array of M complex constellation points, normalized to unit
    average power.
    """
    sqrt_M = int(np.sqrt(M))
    if sqrt_M * sqrt_M != M:
        raise ValueError(f"M must be a perfect square (got {M})")
    # Regular grid from -(sqrt_M-1) to +(sqrt_M-1), step 2
    vals = np.arange(-(sqrt_M - 1), sqrt_M, 2, dtype=float)
    grid = np.array([i + 1j * q for q in vals for i in vals])
    # Normalize to unit average power
    avg_power = np.mean(np.abs(grid) ** 2)
    return grid / np.sqrt(avg_power)


def _qam_mod(bits, M):
    """M-QAM modulation: groups of log2(M) bits -> complex symbols."""
    constellation = _qam_constellation(M)
    k = int(np.log2(M))
    bits = np.asarray(bits, dtype=int)
    # Pad to multiple of k
    pad_len = (-len(bits)) % k
    if pad_len:
        bits = np.append(bits, np.zeros(pad_len, dtype=int))
    num_symbols = len(bits) // k
    symbols = np.empty(num_symbols, dtype=complex)
    for i in range(num_symbols):
        bit_group = bits[i * k:(i + 1) * k]
        idx = int("".join(str(b) for b in bit_group), 2)
        idx = idx % M  # safety clamp
        symbols[i] = constellation[idx]
    return symbols


def _qam_demod(symbols, M):
    """M-QAM demodulation: find nearest constellation point, return bits."""
    constellation = _qam_constellation(M)
    k = int(np.log2(M))
    bits = []
    for s in symbols:
        # Minimum distance detection
        distances = np.abs(s - constellation)
        idx = np.argmin(distances)
        bit_str = format(idx, f"0{k}b")
        bits.extend(int(b) for b in bit_str)
    return np.array(bits, dtype=int)


def simulate_ber_curve(mod_func, demod_func, snr_range, num_bits=100000):
    """Simulate BER across a range of SNR values.

    This is the Monte Carlo approach: generate random bits, modulate them, add
    noise at each SNR level, demodulate, and count errors.  With enough bits
    the result converges to the theoretical curve.

    Parameters
    ----------
    mod_func : callable – modulation function: bits -> complex symbols
    demod_func : callable – demodulation function: noisy symbols -> bits
    snr_range : array-like – SNR values in dB to simulate
    num_bits : int – number of random bits per SNR point (more = smoother)

    Returns
    -------
    ber_values : ndarray – measured BER at each SNR point
    """
    snr_range = np.asarray(snr_range, dtype=float)
    ber_values = np.zeros(len(snr_range))

    tx_bits = np.random.randint(0, 2, num_bits)
    tx_symbols = mod_func(tx_bits)

    for i, snr_db in enumerate(snr_range):
        noisy_symbols, _ = add_awgn(tx_symbols, snr_db)
        rx_bits = demod_func(noisy_symbols)
        # Match lengths (modulator may pad)
        min_len = min(len(tx_bits), len(rx_bits))
        ber, _ = compute_ber(tx_bits[:min_len], rx_bits[:min_len])
        ber_values[i] = ber if ber > 0 else 0.5 / num_bits  # floor for log plot

    return ber_values


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def plot_channel_demo():
    """Generate a multi-panel figure showing channel effects and BER.

    Panels
    ------
    1-4. QPSK constellation at SNR = 30, 15, 5, and 0 dB — shows how noise
         progressively smears the clean constellation points.
    5.   BER vs SNR curves for BPSK, QPSK, 16-QAM, and 64-QAM — both the
         theoretical formulas and Monte Carlo simulations.
    6.   Time domain: a clean BPSK signal vs the same signal corrupted by
         noise, showing how the waveform is distorted.
    7.   Histogram of AWGN noise samples — confirms the Gaussian (bell curve)
         distribution that defines the noise model.
    """
    np.random.seed(42)

    fig = plt.figure(figsize=(16, 14))
    fig.suptitle(
        "Channel Effects & Bit Error Rate (BER)\n"
        "How noise degrades digital signals and causes bit errors",
        fontsize=14,
    )

    # --- Row 1: Constellation degradation (4 panels) -----------------------
    snr_levels = [30, 15, 5, 0]
    num_display_bits = 2000
    display_bits = np.random.randint(0, 2, num_display_bits)
    clean_symbols = _qpsk_mod(display_bits)

    for panel_idx, snr_db in enumerate(snr_levels):
        ax = fig.add_subplot(3, 4, panel_idx + 1)
        noisy_symbols, _ = add_awgn(clean_symbols, snr_db)

        ax.scatter(
            np.real(noisy_symbols), np.imag(noisy_symbols),
            s=4, alpha=0.3, c="tab:blue", edgecolors="none",
        )
        # Mark ideal constellation points
        ideal = _qam_constellation(4)
        ax.scatter(
            np.real(ideal), np.imag(ideal),
            s=80, c="red", marker="x", linewidths=2, zorder=5,
        )
        ax.axhline(0, color="gray", linewidth=0.5)
        ax.axvline(0, color="gray", linewidth=0.5)
        ax.set_xlim(-2.5, 2.5)
        ax.set_ylim(-2.5, 2.5)
        ax.set_aspect("equal")
        ax.set_title(f"QPSK @ SNR = {snr_db} dB", fontsize=10)
        ax.set_xlabel("In-Phase (I)", fontsize=8)
        ax.set_ylabel("Quadrature (Q)", fontsize=8)
        ax.grid(True, alpha=0.3)

        # Compute and display BER for this SNR
        rx_bits = _qpsk_demod(noisy_symbols)
        min_len = min(len(display_bits), len(rx_bits))
        ber, _ = compute_ber(display_bits[:min_len], rx_bits[:min_len])
        quality = "Perfect" if ber == 0 else f"BER = {ber:.4f}"
        ax.text(
            0.05, 0.95, quality, transform=ax.transAxes,
            fontsize=8, verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
        )

    # --- Row 2, left: BER vs SNR curves ------------------------------------
    ax_ber = fig.add_subplot(3, 2, 3)
    snr_range = np.arange(-2, 22, 1.0)
    snr_fine = np.linspace(-2, 22, 200)

    # Theoretical curves
    ax_ber.semilogy(snr_fine, bpsk_theoretical_ber(snr_fine),
                    "b-", linewidth=1.5, label="BPSK (theory)")
    ax_ber.semilogy(snr_fine, qpsk_theoretical_ber(snr_fine),
                    "g--", linewidth=1.5, label="QPSK (theory)")
    ax_ber.semilogy(snr_fine, qam_theoretical_ber(snr_fine, 16),
                    "r-", linewidth=1.5, label="16-QAM (theory)")
    ax_ber.semilogy(snr_fine, qam_theoretical_ber(snr_fine, 64),
                    "m-", linewidth=1.5, label="64-QAM (theory)")

    # Simulated curves (Monte Carlo)
    num_sim_bits = 200000
    ber_bpsk = simulate_ber_curve(_bpsk_mod, _bpsk_demod, snr_range,
                                  num_bits=num_sim_bits)
    ber_qpsk = simulate_ber_curve(_qpsk_mod, _qpsk_demod, snr_range,
                                   num_bits=num_sim_bits)
    ber_16qam = simulate_ber_curve(
        lambda b: _qam_mod(b, 16), lambda s: _qam_demod(s, 16),
        snr_range, num_bits=num_sim_bits,
    )
    ber_64qam = simulate_ber_curve(
        lambda b: _qam_mod(b, 64), lambda s: _qam_demod(s, 64),
        snr_range, num_bits=num_sim_bits,
    )

    ax_ber.semilogy(snr_range, ber_bpsk, "bs", markersize=5, label="BPSK (sim)")
    ax_ber.semilogy(snr_range, ber_qpsk, "g^", markersize=5, label="QPSK (sim)")
    ax_ber.semilogy(snr_range, ber_16qam, "ro", markersize=5, label="16-QAM (sim)")
    ax_ber.semilogy(snr_range, ber_64qam, "md", markersize=5, label="64-QAM (sim)")

    ax_ber.set_xlabel("SNR per bit (dB)")
    ax_ber.set_ylabel("Bit Error Rate (BER)")
    ax_ber.set_title("BER vs SNR — Modulation Comparison")
    ax_ber.legend(fontsize=7, ncol=2, loc="lower left")
    ax_ber.set_ylim(1e-6, 1)
    ax_ber.grid(True, alpha=0.3, which="both")
    ax_ber.text(
        0.98, 0.98,
        "Lower-right = better\n(fewer errors at lower SNR)",
        transform=ax_ber.transAxes, fontsize=7,
        verticalalignment="top", horizontalalignment="right",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.9),
    )

    # --- Row 2, right: Time domain clean vs noisy --------------------------
    ax_time = fig.add_subplot(3, 2, 4)
    num_time_bits = 20
    time_bits = np.random.randint(0, 2, num_time_bits)
    clean_bpsk = _bpsk_mod(time_bits)

    # Upsample for visualization (repeat each symbol for sps samples)
    sps = 50
    clean_wave = np.repeat(clean_bpsk, sps)
    t = np.arange(len(clean_wave))

    noisy_wave, _ = add_awgn(clean_wave, 5.0)

    ax_time.plot(t, clean_wave, "b-", linewidth=1.5, alpha=0.8,
                 label="Clean BPSK signal")
    ax_time.plot(t, noisy_wave, "r-", linewidth=0.8, alpha=0.5,
                 label="With noise (SNR = 5 dB)")
    ax_time.axhline(0, color="gray", linewidth=0.5, linestyle="--")
    ax_time.set_xlabel("Sample index")
    ax_time.set_ylabel("Amplitude")
    ax_time.set_title("Time Domain: Clean vs Noisy Signal")
    ax_time.legend(fontsize=8)
    ax_time.grid(True, alpha=0.3)

    # Annotate bits above the clean waveform
    for i, bit in enumerate(time_bits[:10]):
        ax_time.text(
            i * sps + sps / 2, 1.3, str(bit),
            ha="center", fontsize=7, color="tab:blue",
        )

    # --- Row 3, left: Noise histogram --------------------------------------
    ax_hist = fig.add_subplot(3, 2, 5)
    num_noise_samples = 100000
    noise_samples = np.random.randn(num_noise_samples)

    ax_hist.hist(noise_samples, bins=100, density=True, alpha=0.7,
                 color="tab:blue", edgecolor="none", label="Simulated noise")

    # Overlay theoretical Gaussian curve
    x = np.linspace(-4, 4, 300)
    gaussian_pdf = (1 / np.sqrt(2 * np.pi)) * np.exp(-x ** 2 / 2)
    ax_hist.plot(x, gaussian_pdf, "r-", linewidth=2, label="Gaussian PDF")

    ax_hist.set_xlabel("Noise amplitude")
    ax_hist.set_ylabel("Probability density")
    ax_hist.set_title('AWGN Noise Distribution — "White" = Gaussian Bell Curve')
    ax_hist.legend(fontsize=8)
    ax_hist.grid(True, alpha=0.3)
    ax_hist.text(
        0.98, 0.95,
        "This bell-curve shape is why\nit's called Gaussian noise",
        transform=ax_hist.transAxes, fontsize=7,
        verticalalignment="top", horizontalalignment="right",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.9),
    )

    # --- Row 3, right: Modulation comparison summary -----------------------
    ax_table = fig.add_subplot(3, 2, 6)
    ax_table.axis("off")

    table_data = [
        ["Modulation", "Bits/Symbol", "Spectral\nEfficiency", "Noise\nResilience",
         "Use Case"],
        ["BPSK", "1", "Low", "Best",
         "Deep space, low-SNR links"],
        ["QPSK", "2", "Medium", "Good",
         "Satellite, early LTE"],
        ["16-QAM", "4", "High", "Moderate",
         "Wi-Fi, LTE mid-range"],
        ["64-QAM", "6", "Very High", "Poor",
         "Wi-Fi close-range, cable"],
    ]

    colors = [["#d4e6f1"] * 5]  # header row
    row_colors = ["#ffffff", "#f2f2f2"]
    for i in range(1, 5):
        colors.append([row_colors[i % 2]] * 5)

    table = ax_table.table(
        cellText=table_data, cellColours=colors,
        cellLoc="center", loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.0, 1.6)

    # Bold the header row
    for j in range(5):
        table[0, j].set_text_props(fontweight="bold")

    ax_table.set_title(
        "Modulation Trade-offs: Speed vs Robustness", fontsize=10, pad=15,
    )

    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Interactive widget (for Jupyter notebook)
# ---------------------------------------------------------------------------

def interactive_channel():
    """Launch ipywidgets sliders for exploring channel effects interactively.

    Sliders
    -------
    SNR (dB) : controls how much noise is added. Lower values show the
               constellation smearing and increasing BER.
    Modulation : selects between BPSK, QPSK, 16-QAM, and 64-QAM so you
                 can see how higher-order modulations are more fragile.
    """
    try:
        from ipywidgets import interact, FloatSlider, Dropdown
    except ImportError:
        print("ipywidgets not available — run in Jupyter for interactive mode")
        return

    mod_table = {
        "BPSK": (_bpsk_mod, _bpsk_demod, 2, bpsk_theoretical_ber),
        "QPSK": (_qpsk_mod, _qpsk_demod, 4, qpsk_theoretical_ber),
        "16-QAM": (
            lambda b: _qam_mod(b, 16), lambda s: _qam_demod(s, 16),
            16, lambda snr: qam_theoretical_ber(snr, 16),
        ),
        "64-QAM": (
            lambda b: _qam_mod(b, 64), lambda s: _qam_demod(s, 64),
            64, lambda snr: qam_theoretical_ber(snr, 64),
        ),
    }

    def _update(snr_db=15.0, modulation="QPSK"):
        np.random.seed(0)
        mod_func, demod_func, M, theory_func = mod_table[modulation]

        num_bits = 4000
        tx_bits = np.random.randint(0, 2, num_bits)
        tx_symbols = mod_func(tx_bits)

        noisy_symbols, _ = add_awgn(tx_symbols, snr_db)
        rx_bits = demod_func(noisy_symbols)

        min_len = min(len(tx_bits), len(rx_bits))
        ber, num_errors = compute_ber(tx_bits[:min_len], rx_bits[:min_len])

        fig, axes = plt.subplots(1, 3, figsize=(15, 5))

        # Panel 1: Constellation
        ax = axes[0]
        ax.scatter(
            np.real(noisy_symbols), np.imag(noisy_symbols),
            s=8, alpha=0.3, c="tab:blue", edgecolors="none",
        )
        # Ideal points
        if M == 2:
            ideal = np.array([-1, 1])
            ax.scatter(np.real(ideal), np.imag(ideal), s=80, c="red",
                       marker="x", linewidths=2, zorder=5)
        else:
            ideal = _qam_constellation(M)
            ax.scatter(np.real(ideal), np.imag(ideal), s=80, c="red",
                       marker="x", linewidths=2, zorder=5)
        ax.axhline(0, color="gray", linewidth=0.5)
        ax.axvline(0, color="gray", linewidth=0.5)
        lim = 2.5 if M <= 4 else 1.8
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_title(f"{modulation} Constellation\nSNR = {snr_db:.0f} dB")
        ax.set_xlabel("In-Phase (I)")
        ax.set_ylabel("Quadrature (Q)")
        ax.grid(True, alpha=0.3)

        # Panel 2: BER curve with current operating point
        ax = axes[1]
        snr_fine = np.linspace(-2, 25, 200)
        ax.semilogy(snr_fine, theory_func(snr_fine), "b-", linewidth=1.5,
                     label=f"{modulation} theoretical")
        ax.semilogy(snr_db, max(ber, 1e-6), "ro", markersize=12,
                     label=f"Current: BER = {ber:.2e}", zorder=5)
        ax.axvline(snr_db, color="red", linewidth=0.5, linestyle="--", alpha=0.5)
        ax.set_xlabel("SNR (dB)")
        ax.set_ylabel("BER")
        ax.set_title("BER Curve — You Are Here")
        ax.set_ylim(1e-6, 1)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3, which="both")

        # Panel 3: Error summary text
        ax = axes[2]
        ax.axis("off")
        summary = (
            f"Modulation: {modulation}\n"
            f"SNR: {snr_db:.1f} dB\n"
            f"Bits sent: {min_len:,}\n"
            f"Bit errors: {num_errors:,}\n"
            f"BER: {ber:.2e}\n\n"
        )
        if ber == 0:
            summary += "No errors detected.\nIncrease noise (lower SNR) to see errors."
        elif ber < 1e-3:
            summary += "Acceptable for most applications.\nForward Error Correction can fix these."
        elif ber < 1e-1:
            summary += "Too many errors for reliable data.\nNeed higher SNR or more robust modulation."
        else:
            summary += "Link is essentially broken.\nSignal is buried in noise."
        ax.text(
            0.1, 0.5, summary, transform=ax.transAxes,
            fontsize=11, verticalalignment="center", fontfamily="monospace",
            bbox=dict(boxstyle="round,pad=0.5", facecolor="lightyellow",
                      alpha=0.9),
        )
        ax.set_title("Link Status")

        plt.tight_layout()
        plt.show()

    interact(
        _update,
        snr_db=FloatSlider(
            value=15.0, min=-5.0, max=30.0, step=1.0,
            description="SNR (dB)",
            style={"description_width": "initial"},
        ),
        modulation=Dropdown(
            options=["BPSK", "QPSK", "16-QAM", "64-QAM"],
            value="QPSK",
            description="Modulation",
            style={"description_width": "initial"},
        ),
    )


# ---------------------------------------------------------------------------
# Standalone execution
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    matplotlib.use("Agg")

    fig = plot_channel_demo()
    fig.savefig("channel_ber.png", dpi=150, bbox_inches="tight")
    print("Saved channel_ber.png")

    plt.show()
