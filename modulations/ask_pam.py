"""
ASK (Amplitude Shift Keying) & PAM (Pulse Amplitude Modulation)

These are the simplest ways to send digital data (1s and 0s) using signal
amplitude. Think of them as the "volume knob" approaches to digital
communication.

ASK — Amplitude Shift Keying:
  The most intuitive digital modulation: to send a 1, transmit a carrier wave;
  to send a 0, turn it off (or make it quieter). This on/off variant is called
  OOK (On-Off Keying) and it's essentially how a flashlight sends Morse code —
  light on = 1, light off = 0. More generally, ASK maps each digital symbol to
  a different carrier amplitude level.

PAM — Pulse Amplitude Modulation:
  Instead of switching a carrier on/off, PAM works directly at baseband (no
  carrier wave). Each digital symbol becomes a pulse whose *height* encodes
  the data. A 2-PAM signal has two heights (+1 and -1) — one bit per symbol.
  A 4-PAM signal has four heights (-3, -1, +1, +3) — two bits per symbol.
  An 8-PAM signal has eight heights — three bits per symbol, and so on.

  More levels (higher M) means you pack more bits into each symbol period, so
  the data rate goes up. The catch: the amplitude levels get closer together,
  making it harder for the receiver to tell them apart in noise. This is the
  fundamental throughput-vs-reliability tradeoff in digital communications.

Eye diagram:
  A diagnostic tool that overlays many consecutive symbol periods on top of
  each other. If the "eye" is wide open, the receiver has plenty of margin to
  distinguish symbol levels. As noise or distortion increases, the eye closes
  — and bit errors start piling up. Engineers use eye diagrams the way a
  doctor uses an EKG: a quick visual health check of the signal.

Why this matters:
  ASK and PAM are building blocks for more complex modulation schemes. QAM
  (Quadrature Amplitude Modulation), used in Wi-Fi, LTE, and cable modems,
  is essentially two PAM signals combined on carriers 90 degrees apart.
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Core functions (importable by notebook)
# ---------------------------------------------------------------------------

def ask_modulate(bits, carrier_freq, fs, sps):
    """Generate an OOK (On-Off Keying) ASK signal from a bit sequence.

    Each bit is held for *sps* samples. A '1' bit transmits a carrier sine
    wave at full amplitude; a '0' bit transmits silence.

    Parameters
    ----------
    bits : array-like – sequence of 0s and 1s
    carrier_freq : float – carrier frequency in Hz
    fs : float – sample rate in samples/second
    sps : int – samples per symbol (controls symbol duration)

    Returns
    -------
    t : ndarray – time vector (seconds)
    signal : ndarray – ASK-modulated waveform
    message : ndarray – rectangular pulse train of the original bits
    """
    bits = np.asarray(bits, dtype=float)
    num_samples = len(bits) * sps
    t = np.arange(num_samples) / fs

    # Rectangular pulse train: repeat each bit value for sps samples
    message = np.repeat(bits, sps)

    # Carrier wave
    carrier = np.sin(2 * np.pi * carrier_freq * t)

    # ASK = carrier * message (OOK: carrier on when bit=1, off when bit=0)
    signal = carrier * message
    return t, signal, message


def pam_modulate(bits, M, sps):
    """Generate an M-ary PAM baseband signal.

    Groups the input bits into symbols of log2(M) bits each, maps each symbol
    to one of M equally-spaced amplitude levels, and produces a pulse train.

    Parameters
    ----------
    bits : array-like – sequence of 0s and 1s (length must be divisible
           by log2(M))
    M : int – constellation size (number of amplitude levels, must be a
        power of 2)
    sps : int – samples per symbol

    Returns
    -------
    signal : ndarray – baseband PAM waveform
    symbols : ndarray – integer symbol values before amplitude mapping
    levels : ndarray – the M amplitude levels used
    """
    bits = np.asarray(bits, dtype=int)
    k = int(np.log2(M))  # bits per symbol
    if len(bits) % k != 0:
        raise ValueError(
            f"Bit length ({len(bits)}) must be divisible by log2(M)={k}"
        )

    # Reshape bits into groups of k and convert to symbol indices
    bit_groups = bits.reshape(-1, k)
    symbols = np.zeros(len(bit_groups), dtype=int)
    for i, group in enumerate(bit_groups):
        # MSB-first binary to integer
        symbols[i] = int("".join(str(b) for b in group), 2)

    # Map symbol indices to equally-spaced amplitude levels
    # For M=4: levels = [-3, -1, +1, +3]
    levels = np.arange(M) * 2 - (M - 1)
    amplitudes = levels[symbols].astype(float)

    # Upsample: hold each amplitude for sps samples
    signal = np.repeat(amplitudes, sps)
    return signal, symbols, levels


def generate_eye_diagram(signal, sps, num_traces=50, offset=0):
    """Extract overlapping symbol-period segments for an eye diagram.

    An eye diagram overlays consecutive two-symbol-wide windows of the signal.
    Where the traces cluster tightly, the signal is clean; where they spread
    apart, noise or distortion is degrading quality.

    Parameters
    ----------
    signal : ndarray – the waveform to analyze
    sps : int – samples per symbol
    num_traces : int – number of overlaid segments to extract
    offset : int – starting sample offset

    Returns
    -------
    traces : list of ndarray – each element is a 2-symbol-wide segment
    """
    window = 2 * sps  # two symbol periods per trace
    traces = []
    start = offset
    for _ in range(num_traces):
        if start + window > len(signal):
            break
        traces.append(signal[start:start + window])
        start += sps  # advance by one symbol period
    return traces


def add_awgn(signal, snr_db):
    """Add white Gaussian noise to a signal at a specified SNR.

    AWGN (Additive White Gaussian Noise) is the standard noise model in
    communications — it represents the thermal noise present in all
    electronic receivers.

    Parameters
    ----------
    signal : ndarray – clean input signal
    snr_db : float – desired signal-to-noise ratio in decibels

    Returns
    -------
    noisy : ndarray – signal + noise
    """
    sig_power = np.mean(signal ** 2)
    snr_linear = 10 ** (snr_db / 10)
    noise_power = sig_power / snr_linear
    noise = np.sqrt(noise_power) * np.random.randn(len(signal))
    return signal + noise


def compute_ber(bits_tx, bits_rx):
    """Compute Bit Error Rate between transmitted and received bit sequences."""
    bits_tx = np.asarray(bits_tx)
    bits_rx = np.asarray(bits_rx)
    return np.mean(bits_tx != bits_rx)


def pam_demodulate(signal, M, sps):
    """Demodulate a noisy M-PAM signal by sampling at symbol centers and
    mapping to the nearest valid amplitude level.

    Parameters
    ----------
    signal : ndarray – received (possibly noisy) PAM waveform
    M : int – constellation size
    sps : int – samples per symbol

    Returns
    -------
    bits_rx : ndarray – recovered bit sequence
    """
    k = int(np.log2(M))
    levels = np.arange(M) * 2 - (M - 1)

    # Sample at the center of each symbol period
    sample_points = np.arange(sps // 2, len(signal), sps)
    samples = signal[sample_points]

    # Find nearest constellation level for each sample
    symbol_indices = np.array([
        np.argmin(np.abs(levels - s)) for s in samples
    ])

    # Convert symbol indices back to bits (MSB-first)
    bits_rx = []
    for idx in symbol_indices:
        bit_str = format(idx, f"0{k}b")
        bits_rx.extend(int(b) for b in bit_str)
    return np.array(bits_rx)


def simulate_ber_vs_snr(M, snr_range, num_bits=30000, sps=8):
    """Run a Monte Carlo BER simulation for M-PAM over a range of SNR values.

    Parameters
    ----------
    M : int – PAM constellation size
    snr_range : array-like – SNR values in dB to simulate
    num_bits : int – number of random bits per trial
    sps : int – samples per symbol

    Returns
    -------
    ber_values : ndarray – measured BER at each SNR point
    """
    k = int(np.log2(M))
    # Ensure num_bits is divisible by k
    num_bits = (num_bits // k) * k
    bits_tx = np.random.randint(0, 2, num_bits)

    signal_clean, _, _ = pam_modulate(bits_tx, M, sps)

    ber_values = np.zeros(len(snr_range))
    for i, snr_db in enumerate(snr_range):
        noisy = add_awgn(signal_clean, snr_db)
        bits_rx = pam_demodulate(noisy, M, sps)
        ber_values[i] = compute_ber(bits_tx, bits_rx)

    return ber_values


# ---------------------------------------------------------------------------
# Plotting — static demo figure
# ---------------------------------------------------------------------------

def plot_ask_pam_demo():
    """Generate a multi-panel figure demonstrating ASK and PAM concepts.

    Panels:
      1. Digital message (bit sequence)
      2. ASK modulated signal (on-off keying)
      3. 4-PAM signal with symbol boundaries
      4. Eye diagram for 4-PAM with noise
      5. BER vs SNR for different PAM orders
    """
    np.random.seed(42)

    # -- Shared parameters --
    bits = np.array([1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1])
    sps = 100  # samples per symbol (smooth waveforms for plotting)
    carrier_freq = 50.0  # Hz
    fs = 5000.0  # sample rate

    fig, axes = plt.subplots(5, 1, figsize=(12, 14))
    fig.suptitle(
        "ASK & PAM  —  Digital Modulation via Amplitude",
        fontsize=14, fontweight="bold",
    )

    # --- Panel 1: Digital message ---
    ax = axes[0]
    t_bits = np.arange(len(bits) * sps) / fs
    bit_waveform = np.repeat(bits.astype(float), sps)
    ax.plot(t_bits * 1000, bit_waveform, color="tab:blue", linewidth=1.5)
    ax.set_ylim(-0.3, 1.5)
    ax.set_title("Digital Message (Bit Sequence)")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Bit Value")
    # Mark bit values at center of each symbol
    for i, b in enumerate(bits):
        t_center = (i + 0.5) * sps / fs * 1000
        ax.text(t_center, 1.2, str(b), ha="center", fontsize=9,
                fontweight="bold", color="tab:blue")
    ax.grid(True, alpha=0.3)

    # --- Panel 2: ASK (OOK) modulated signal ---
    ax = axes[1]
    t_ask, ask_signal, _ = ask_modulate(bits, carrier_freq, fs, sps)
    ax.plot(t_ask * 1000, ask_signal, color="tab:orange", linewidth=0.8)
    ax.set_title("ASK Modulated Signal (On-Off Keying)")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude")
    ax.grid(True, alpha=0.3)
    # Add shading for '0' bits
    for i, b in enumerate(bits):
        t_start = i * sps / fs * 1000
        t_end = (i + 1) * sps / fs * 1000
        if b == 0:
            ax.axvspan(t_start, t_end, alpha=0.1, color="red")

    # --- Panel 3: 4-PAM with symbol boundaries ---
    ax = axes[2]
    M = 4
    pam_bits = np.array([1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1])  # 12 bits -> 6 symbols
    pam_signal, pam_symbols, pam_levels = pam_modulate(pam_bits, M, sps)
    k = int(np.log2(M))
    t_pam = np.arange(len(pam_signal)) / fs
    ax.plot(t_pam * 1000, pam_signal, color="tab:green", linewidth=1.5)
    ax.set_title(f"{M}-PAM Signal (2 bits per symbol)")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude Level")
    # Draw horizontal lines at each amplitude level
    for level in pam_levels:
        ax.axhline(y=level, color="gray", linestyle="--", alpha=0.4)
    # Mark symbol boundaries and label bit groups
    num_symbols = len(pam_bits) // k
    for i in range(num_symbols):
        t_boundary = i * sps / fs * 1000
        ax.axvline(x=t_boundary, color="gray", linestyle=":", alpha=0.5)
        # Label the bit group above the signal
        bit_label = "".join(str(b) for b in pam_bits[i * k:(i + 1) * k])
        t_center = (i + 0.5) * sps / fs * 1000
        ax.text(t_center, pam_levels[-1] + 0.8, bit_label, ha="center",
                fontsize=8, color="tab:green", fontweight="bold")
    ax.grid(True, alpha=0.3)

    # --- Panel 4: Eye diagram for noisy 4-PAM ---
    ax = axes[3]
    # Generate a longer signal for a meaningful eye diagram
    long_bits = np.random.randint(0, 2, 400)
    sps_eye = 32
    pam_long, _, _ = pam_modulate(long_bits, M, sps_eye)
    pam_noisy = add_awgn(pam_long, snr_db=15)
    traces = generate_eye_diagram(pam_noisy, sps_eye, num_traces=80)
    t_eye = np.linspace(0, 2, 2 * sps_eye)
    for trace in traces:
        ax.plot(t_eye, trace, color="tab:purple", alpha=0.15, linewidth=0.8)
    ax.set_title("Eye Diagram  (4-PAM, SNR = 15 dB)")
    ax.set_xlabel("Time (symbol periods)")
    ax.set_ylabel("Amplitude")
    ax.axvline(x=1.0, color="red", linestyle="--", alpha=0.5, label="Optimal sample point")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # --- Panel 5: BER vs SNR for different PAM orders ---
    ax = axes[4]
    snr_range = np.arange(0, 26, 1)
    for M_order, color, marker in [(2, "tab:blue", "o"),
                                    (4, "tab:green", "s"),
                                    (8, "tab:red", "^")]:
        ber = simulate_ber_vs_snr(M_order, snr_range, num_bits=30000, sps=8)
        # Replace zeros with a small value for log plotting
        ber_plot = np.where(ber > 0, ber, 1e-6)
        ax.semilogy(snr_range, ber_plot, marker=marker, markersize=4,
                    linewidth=1.2, color=color,
                    label=f"{M_order}-PAM ({int(np.log2(M_order))} bit/sym)")
    ax.set_title("BER vs SNR  (more levels = more errors at same SNR)")
    ax.set_xlabel("SNR (dB)")
    ax.set_ylabel("Bit Error Rate")
    ax.set_ylim(1e-5, 1)
    ax.legend(fontsize=9)
    ax.grid(True, which="both", alpha=0.3)

    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Interactive widget (for Jupyter notebook)
# ---------------------------------------------------------------------------

def interactive_ask_pam():
    """Launch ipywidgets sliders for exploring ASK and PAM interactively.

    Sliders:
      - SNR (dB): controls noise level added to the PAM signal
      - PAM order (M): selects 2-PAM, 4-PAM, or 8-PAM
    """
    try:
        from ipywidgets import interact, IntSlider, Dropdown
    except ImportError:
        print("ipywidgets not available -- run in Jupyter for interactive mode")
        return

    def _update(snr_db=15, M=4):
        np.random.seed(0)
        k = int(np.log2(M))
        num_bits = 400
        num_bits = (num_bits // k) * k
        sps = 32
        bits = np.random.randint(0, 2, num_bits)

        # Generate PAM signal and add noise
        signal_clean, symbols, levels = pam_modulate(bits, M, sps)
        signal_noisy = add_awgn(signal_clean, snr_db)

        # Demodulate and compute BER
        bits_rx = pam_demodulate(signal_noisy, M, sps)
        ber = compute_ber(bits, bits_rx)

        fig, axes = plt.subplots(1, 3, figsize=(14, 4))
        fig.suptitle(
            f"{M}-PAM  |  SNR = {snr_db} dB  |  BER = {ber:.4f}",
            fontsize=12, fontweight="bold",
        )

        # Clean vs noisy signal (first 5 symbols)
        ax = axes[0]
        show_samples = 5 * sps
        t = np.arange(show_samples) / (sps * 1.0)
        ax.plot(t, signal_clean[:show_samples], linewidth=2,
                label="Clean", color="tab:blue")
        ax.plot(t, signal_noisy[:show_samples], linewidth=0.8,
                alpha=0.7, label="Noisy", color="tab:orange")
        ax.set_title("PAM Waveform")
        ax.set_xlabel("Symbol periods")
        ax.set_ylabel("Amplitude")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        # Eye diagram
        ax = axes[1]
        traces = generate_eye_diagram(signal_noisy, sps, num_traces=60)
        t_eye = np.linspace(0, 2, 2 * sps)
        for trace in traces:
            ax.plot(t_eye, trace, color="tab:purple", alpha=0.15,
                    linewidth=0.8)
        ax.set_title("Eye Diagram")
        ax.set_xlabel("Symbol periods")
        ax.set_ylabel("Amplitude")
        ax.grid(True, alpha=0.3)

        # Constellation (received sample points)
        ax = axes[2]
        sample_points = np.arange(sps // 2, len(signal_noisy), sps)
        rx_samples = signal_noisy[sample_points]
        ax.hist(rx_samples, bins=50, orientation="horizontal",
                color="tab:green", alpha=0.7, edgecolor="black",
                linewidth=0.3)
        for level in levels:
            ax.axhline(y=level, color="red", linestyle="--", alpha=0.6)
        ax.set_title("Received Amplitude Histogram")
        ax.set_xlabel("Count")
        ax.set_ylabel("Amplitude")
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

    interact(
        _update,
        snr_db=IntSlider(value=15, min=0, max=30, step=1,
                         description="SNR (dB)"),
        M=Dropdown(options=[2, 4, 8], value=4, description="PAM Order (M)"),
    )


# ---------------------------------------------------------------------------
# Standalone execution
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    matplotlib.use("Agg")
    fig = plot_ask_pam_demo()
    fig.savefig("ask_pam.png", dpi=150)
    print("Saved ask_pam.png")

    plt.show()
