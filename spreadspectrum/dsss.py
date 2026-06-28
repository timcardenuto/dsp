"""
Direct Sequence Spread Spectrum (DSSS)

Imagine you want to send a message in a crowded room. Instead of shouting one
clear word (which anyone can hear and a jammer can drown out), you encode your
message into a long, rapid, noise-like pattern that only your receiver knows
how to decode. That is the core idea behind DSSS.

How it works:
  1. SPREAD — Take your original data signal and multiply every bit by a much
     faster pseudo-random (PN) sequence. Each data bit gets replaced by many
     fast "chips," spreading the signal's energy across a much wider bandwidth.
     In the frequency domain the tall, narrow peak becomes a low, wide hump
     that can even sit below the noise floor.

  2. TRANSMIT — The spread signal looks like noise to anyone who doesn't know
     the PN sequence. It is hard to detect, hard to jam, and multiple users
     can share the same band with different PN codes (this is how CDMA works).

  3. DESPREAD — The receiver multiplies the received signal by the exact same
     PN sequence. Because the PN sequence XORed with itself equals all-ones
     (i.e., it cancels out), the original data pops back out. Any interfering
     signal that does NOT match the PN code stays spread and appears as
     low-level noise.

Key concept — Processing Gain:
  Processing gain = spread bandwidth / original bandwidth
                  = number of PN chips per data bit
  A processing gain of 10 means the system can tolerate 10x more interference
  power than a conventional system. Higher gain = more jam resistance, but you
  need more bandwidth.

Why use DSSS?
  - Jam resistance: spreading makes it very hard to block the signal
  - Low probability of intercept: signal hides below the noise floor
  - Multiple access (CDMA): many users share one band with unique PN codes
  - Multipath resistance: spreading reduces sensitivity to reflections
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Core functions (importable by notebook)
# ---------------------------------------------------------------------------

def generate_pn_sequence(length, seed=None):
    """Generate a pseudo-random noise (PN) sequence of +1 and -1 values.

    A PN sequence is a deterministic but noise-like binary pattern. In real
    systems these come from Linear Feedback Shift Registers (LFSRs); here we
    use a seeded random generator for simplicity.

    Parameters
    ----------
    length : int – number of chips (fast bits) in the sequence
    seed : int or None – random seed for reproducibility

    Returns
    -------
    pn : ndarray of {+1, -1} with shape (length,)
    """
    rng = np.random.default_rng(seed)
    # Generate 0/1 then map to +1/-1 (bipolar representation)
    bits = rng.integers(0, 2, size=length)
    pn = 2 * bits - 1  # 0 -> -1, 1 -> +1
    return pn


def dsss_spread(data, pn_sequence):
    """Spread a data signal using a PN sequence.

    Each data sample is multiplied by every chip in the PN sequence, producing
    an output that is len(data) * len(pn_sequence) samples long.

    Think of it like replacing every single letter of your message with a whole
    paragraph of coded text — the message gets much longer but much harder to
    read without the codebook.

    Parameters
    ----------
    data : array-like of {+1, -1} – bipolar data symbols
    pn_sequence : array-like of {+1, -1} – PN spreading code

    Returns
    -------
    spread : ndarray – spread signal, length = len(data) * len(pn_sequence)
    """
    data = np.asarray(data, dtype=float)
    pn = np.asarray(pn_sequence, dtype=float)
    chips_per_bit = len(pn)

    # Each data bit multiplied by the full PN sequence
    spread = np.repeat(data, chips_per_bit) * np.tile(pn, len(data))
    return spread


def dsss_despread(signal, pn_sequence):
    """Despread a received signal using the same PN sequence.

    Multiplies the received signal chip-by-chip with the PN sequence, then
    integrates (sums) over each group of chips to recover one data symbol.
    This "integrate and dump" step is what collapses the spread energy back
    into a single strong value per bit.

    Parameters
    ----------
    signal : array-like – received spread signal
    pn_sequence : array-like of {+1, -1} – PN code (must match transmitter)

    Returns
    -------
    recovered : ndarray – recovered bipolar data symbols (+1 or -1)
    """
    signal = np.asarray(signal, dtype=float)
    pn = np.asarray(pn_sequence, dtype=float)
    chips_per_bit = len(pn)
    num_bits = len(signal) // chips_per_bit

    # Multiply by PN sequence then sum each bit-length block
    despread = signal[:num_bits * chips_per_bit] * np.tile(pn, num_bits)
    # Reshape into (num_bits, chips_per_bit) and sum each row
    sums = despread.reshape(num_bits, chips_per_bit).sum(axis=1)
    # Decision: positive sum -> +1, negative sum -> -1
    recovered = np.sign(sums)
    return recovered


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def _compute_psd(signal, fs):
    """Compute power spectral density in dB (helper for plots)."""
    N = len(signal)
    Y = np.fft.fft(signal)
    psd = (1 / (fs * N)) * np.abs(Y[:N // 2]) ** 2
    psd[1:] *= 2
    freqs = np.fft.fftfreq(N, 1 / fs)[:N // 2]
    psd_db = 10 * np.log10(psd + 1e-12)
    return freqs, psd_db


def plot_dsss_demo(num_bits=10, chips_per_bit=16, seed=42):
    """Generate the multi-panel DSSS overview figure.

    Shows the complete spread-spectrum pipeline: original data, PN sequence,
    spreading, spectrum comparison, and despreading recovery.

    Parameters
    ----------
    num_bits : int – number of data bits to transmit
    chips_per_bit : int – PN sequence length (= processing gain)
    seed : int – random seed for reproducibility
    """
    rng = np.random.default_rng(seed)

    # --- Generate data and PN sequence ---
    data_bits = 2 * rng.integers(0, 2, size=num_bits) - 1  # bipolar +1/-1
    pn = generate_pn_sequence(chips_per_bit, seed=seed + 1)

    # --- Spread ---
    spread_signal = dsss_spread(data_bits, pn)

    # --- Despread (recover) ---
    recovered = dsss_despread(spread_signal, pn)

    # --- Build time axes ---
    # For data: each bit occupies `chips_per_bit` chip-periods
    data_upsampled = np.repeat(data_bits, chips_per_bit)
    chip_indices = np.arange(len(spread_signal))
    pn_tiled = np.tile(pn, num_bits)

    # --- Spectrum comparison ---
    # Use a simple sample rate of 1 sample per chip
    fs_chip = 1.0
    # Original narrowband data (upsampled to same length for fair comparison)
    freqs_data, psd_data = _compute_psd(data_upsampled, fs_chip)
    freqs_spread, psd_spread = _compute_psd(spread_signal, fs_chip)

    # --- Plot ---
    fig, axes = plt.subplots(3, 2, figsize=(14, 10))
    fig.suptitle("Direct Sequence Spread Spectrum (DSSS)", fontsize=14)

    # Panel 1: Original data signal
    ax = axes[0, 0]
    ax.step(chip_indices, data_upsampled, where="post", color="tab:blue",
            linewidth=1.5)
    ax.set_title("Original Data Signal")
    ax.set_xlabel("Chip index")
    ax.set_ylabel("Amplitude")
    ax.set_ylim(-1.5, 1.5)
    ax.set_yticks([-1, 0, 1])
    ax.grid(True, alpha=0.3)
    # Annotate bit values
    for i, bit in enumerate(data_bits):
        ax.text((i + 0.5) * chips_per_bit, 1.2,
                f"{'+1' if bit > 0 else '-1'}", ha="center", fontsize=7,
                color="tab:blue")

    # Panel 2: PN spreading sequence (one period, then tiled)
    ax = axes[0, 1]
    ax.step(chip_indices, pn_tiled, where="post", color="tab:orange",
            linewidth=0.8)
    ax.set_title("PN Spreading Sequence (repeats each bit)")
    ax.set_xlabel("Chip index")
    ax.set_ylabel("Amplitude")
    ax.set_ylim(-1.5, 1.5)
    ax.set_yticks([-1, 0, 1])
    ax.grid(True, alpha=0.3)

    # Panel 3: Spread signal (data * PN)
    ax = axes[1, 0]
    ax.step(chip_indices, spread_signal, where="post", color="tab:green",
            linewidth=0.8)
    ax.set_title("Spread Signal (Data x PN)")
    ax.set_xlabel("Chip index")
    ax.set_ylabel("Amplitude")
    ax.set_ylim(-1.5, 1.5)
    ax.set_yticks([-1, 0, 1])
    ax.grid(True, alpha=0.3)

    # Panel 4: Spectrum comparison — before and after spreading
    ax = axes[1, 1]
    ax.plot(freqs_data, psd_data, color="tab:blue", alpha=0.8,
            label="Original (narrowband)")
    ax.plot(freqs_spread, psd_spread, color="tab:green", alpha=0.8,
            label="Spread (wideband)")
    ax.set_title("Power Spectrum: Before vs After Spreading")
    ax.set_xlabel("Normalized Frequency")
    ax.set_ylabel("Power (dB)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel 5: Recovered / despread signal
    ax = axes[2, 0]
    recovered_upsampled = np.repeat(recovered, chips_per_bit)
    ax.step(chip_indices, recovered_upsampled, where="post",
            color="tab:purple", linewidth=1.5)
    ax.set_title("Despread / Recovered Signal")
    ax.set_xlabel("Chip index")
    ax.set_ylabel("Amplitude")
    ax.set_ylim(-1.5, 1.5)
    ax.set_yticks([-1, 0, 1])
    ax.grid(True, alpha=0.3)
    # Show match status
    match = np.array_equal(data_bits, recovered)
    ax.text(0.98, 0.95,
            f"{'All bits recovered correctly' if match else 'ERRORS in recovery'}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9,
            color="green" if match else "red",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    # Panel 6: Despreading with WRONG PN code (demonstrates rejection)
    ax = axes[2, 1]
    wrong_pn = generate_pn_sequence(chips_per_bit, seed=seed + 99)
    wrong_recovered = dsss_despread(spread_signal, wrong_pn)
    wrong_upsampled = np.repeat(wrong_recovered, chips_per_bit)
    ax.step(chip_indices, wrong_upsampled, where="post", color="tab:red",
            linewidth=1.5)
    ax.set_title("Despread with WRONG PN Code (signal rejected)")
    ax.set_xlabel("Chip index")
    ax.set_ylabel("Amplitude")
    ax.set_ylim(-1.5, 1.5)
    ax.set_yticks([-1, 0, 1])
    ax.grid(True, alpha=0.3)
    wrong_match = np.array_equal(data_bits, wrong_recovered)
    ax.text(0.98, 0.95,
            f"{'Match (unlikely!)' if wrong_match else 'Data NOT recovered (as expected)'}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9,
            color="green" if wrong_match else "red",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Interactive widget (for Jupyter notebook)
# ---------------------------------------------------------------------------

def interactive_dsss():
    """Launch ipywidgets sliders for exploring DSSS behavior.

    Sliders
    -------
    - chips_per_bit : PN sequence length (= processing gain). Higher values
      spread the signal more, making it harder to jam but using more bandwidth.
    - snr_db : Signal-to-noise ratio. Lower values simulate a noisier channel.
      With enough processing gain, DSSS can recover data even at very low SNR.
    """
    try:
        from ipywidgets import interact, IntSlider, FloatSlider
    except ImportError:
        print("ipywidgets not available -- run in Jupyter for interactive mode")
        return

    def _update(chips_per_bit=16, snr_db=10.0):
        num_bits = 10
        seed = 42
        rng = np.random.default_rng(seed)

        # Generate data and PN sequence
        data_bits = 2 * rng.integers(0, 2, size=num_bits) - 1
        pn = generate_pn_sequence(chips_per_bit, seed=seed + 1)

        # Spread
        spread_signal = dsss_spread(data_bits, pn)

        # Add noise
        signal_power = np.mean(spread_signal ** 2)
        snr_linear = 10 ** (snr_db / 10)
        noise_power = signal_power / snr_linear
        noise = np.sqrt(noise_power) * np.random.default_rng(0).standard_normal(
            len(spread_signal))
        received = spread_signal + noise

        # Despread noisy signal
        recovered = dsss_despread(received, pn)

        # Bit error rate
        errors = np.sum(data_bits != recovered)
        ber = errors / num_bits

        # Build time axes
        chip_indices = np.arange(len(spread_signal))
        data_upsampled = np.repeat(data_bits, chips_per_bit)
        recovered_upsampled = np.repeat(recovered, chips_per_bit)

        # Spectra
        fs_chip = 1.0
        freqs_data, psd_data = _compute_psd(data_upsampled, fs_chip)
        freqs_spread, psd_spread = _compute_psd(received, fs_chip)

        fig, axes = plt.subplots(2, 2, figsize=(14, 8))
        fig.suptitle(
            f"DSSS Interactive  |  Processing Gain = {chips_per_bit}  |  "
            f"SNR = {snr_db:.0f} dB  |  BER = {ber:.0%}",
            fontsize=13)

        # Original data
        ax = axes[0, 0]
        ax.step(chip_indices, data_upsampled, where="post", color="tab:blue",
                linewidth=1.5)
        ax.set_title("Original Data")
        ax.set_ylabel("Amplitude")
        ax.set_ylim(-1.5, 1.5)
        ax.grid(True, alpha=0.3)

        # Received (spread + noise)
        ax = axes[0, 1]
        ax.plot(chip_indices, received, color="tab:green", linewidth=0.5,
                alpha=0.8)
        ax.set_title("Received Signal (spread + noise)")
        ax.set_ylabel("Amplitude")
        ax.grid(True, alpha=0.3)

        # Spectrum comparison
        ax = axes[1, 0]
        ax.plot(freqs_data, psd_data, color="tab:blue", alpha=0.8,
                label="Original")
        ax.plot(freqs_spread, psd_spread, color="tab:green", alpha=0.8,
                label="Received (spread)")
        ax.set_title("Power Spectrum Comparison")
        ax.set_xlabel("Normalized Frequency")
        ax.set_ylabel("Power (dB)")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        # Recovered data
        ax = axes[1, 1]
        ax.step(chip_indices, recovered_upsampled, where="post",
                color="tab:purple", linewidth=1.5)
        ax.set_title(f"Recovered Data  ({errors} bit errors out of {num_bits})")
        ax.set_xlabel("Chip index")
        ax.set_ylabel("Amplitude")
        ax.set_ylim(-1.5, 1.5)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

    interact(
        _update,
        chips_per_bit=IntSlider(value=16, min=4, max=64, step=4,
                                description="Chips/bit (PG)"),
        snr_db=FloatSlider(value=10.0, min=-10.0, max=30.0, step=1.0,
                           description="SNR (dB)"),
    )


# ---------------------------------------------------------------------------
# Standalone execution
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    matplotlib.use("Agg")
    fig = plot_dsss_demo()
    fig.savefig("dsss_demo.png", dpi=150)
    print("Saved dsss_demo.png")

    plt.show()
