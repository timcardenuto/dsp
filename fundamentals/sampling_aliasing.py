"""
Sampling & Aliasing

In the real world, signals are continuous — a microphone produces a smoothly
varying voltage, not a list of numbers. To process signals on a computer, we
must "sample" them: measure the signal's value at evenly spaced points in time
and store each measurement as a number. The rate at which we take these
measurements is called the sample rate (fs), measured in samples per second (Hz).

The Nyquist theorem tells us the minimum sample rate needed to faithfully
capture a signal: fs must be at least TWICE the highest frequency present in
the signal (fs >= 2 * fmax). This critical threshold (fs/2) is called the
Nyquist frequency.

What happens if we sample too slowly? The answer is aliasing — the signal
"folds" and appears at a completely wrong frequency. Imagine a car wheel in a
movie: if the camera frame rate is too low, the wheel appears to spin backward.
That's aliasing. A 900 Hz tone sampled at 1000 Hz looks identical to a 100 Hz
tone — the high frequency masquerades as a low one.

Reconstruction is the reverse process: given a set of discrete samples, we
rebuild an approximation of the original continuous signal. When the Nyquist
criterion is satisfied, perfect reconstruction is theoretically possible using
sinc interpolation (each sample contributes a sinc-shaped pulse, and their sum
recreates the continuous waveform).
"""

import numpy as np
from scipy import signal as sp_signal
import matplotlib
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Core functions (importable by notebook)
# ---------------------------------------------------------------------------

def sample_signal(freq, duration, fs, amplitude=1.0, phase=0.0):
    """Sample a continuous sinusoid at a given sample rate.

    Parameters
    ----------
    freq : float – frequency of the sinusoid in Hz
    duration : float – length in seconds
    fs : float – sample rate in samples/second
    amplitude : float – peak amplitude
    phase : float – phase offset in radians

    Returns
    -------
    t_continuous : ndarray – dense time vector (simulating continuous)
    x_continuous : ndarray – dense signal values
    t_sampled : ndarray – sample time instants
    x_sampled : ndarray – sampled values
    """
    # "Continuous" signal — use a very high sample rate to approximate
    fs_continuous = max(50 * freq, 10000)
    t_continuous = np.arange(0, duration, 1 / fs_continuous)
    x_continuous = amplitude * np.sin(2 * np.pi * freq * t_continuous + phase)

    # Sampled version
    t_sampled = np.arange(0, duration, 1 / fs)
    x_sampled = amplitude * np.sin(2 * np.pi * freq * t_sampled + phase)

    return t_continuous, x_continuous, t_sampled, x_sampled


def demonstrate_aliasing(freq, fs):
    """Show what frequency a sampled signal actually appears at.

    When fs < 2*freq, the signal aliases to a different frequency.
    The apparent (aliased) frequency is: |freq - round(freq/fs) * fs|

    Parameters
    ----------
    freq : float – true signal frequency in Hz
    fs : float – sample rate in Hz

    Returns
    -------
    aliased_freq : float – the frequency the signal appears at after sampling
    is_aliased : bool – True if the signal is aliased (fs < 2*freq)
    """
    nyquist = fs / 2.0
    is_aliased = freq > nyquist

    # Compute aliased frequency by folding into [0, fs/2]
    aliased_freq = freq % fs
    if aliased_freq > nyquist:
        aliased_freq = fs - aliased_freq

    return aliased_freq, is_aliased


def reconstruct_signal(t_sampled, x_sampled, fs, t_recon):
    """Reconstruct a continuous signal from samples using sinc interpolation.

    This implements the Whittaker-Shannon interpolation formula:
      x(t) = sum_n [ x[n] * sinc((t - n/fs) * fs) ]

    When the Nyquist criterion was met during sampling, this gives perfect
    reconstruction. When it was violated, you get the aliased signal back.

    Parameters
    ----------
    t_sampled : ndarray – sample time instants
    x_sampled : ndarray – sampled values
    fs : float – sample rate used during sampling
    t_recon : ndarray – time points at which to reconstruct

    Returns
    -------
    x_recon : ndarray – reconstructed signal values
    """
    # Sinc interpolation: each sample contributes a sinc pulse
    T = 1.0 / fs
    x_recon = np.zeros_like(t_recon, dtype=float)
    for n, (tn, xn) in enumerate(zip(t_sampled, x_sampled)):
        x_recon += xn * np.sinc((t_recon - tn) / T)
    return x_recon


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def plot_sampling_demo(freq=5.0, fs_good=40.0, fs_bad=8.0, duration=1.0):
    """Generate a four-panel figure demonstrating sampling and aliasing.

    Parameters
    ----------
    freq : float – signal frequency in Hz (default 5 Hz for clear visuals)
    fs_good : float – a sample rate above Nyquist (default 40 Hz)
    fs_bad : float – a sample rate below Nyquist (default 8 Hz)
    duration : float – signal duration in seconds
    """
    # --- Generate signals ---
    t_cont, x_cont, t_good, x_good = sample_signal(freq, duration, fs_good)
    _, _, t_bad, x_bad = sample_signal(freq, duration, fs_bad)

    # Aliased frequency
    aliased_freq, _ = demonstrate_aliasing(freq, fs_bad)

    # Reconstruct from undersampled data
    t_recon = t_cont  # reconstruct onto dense time grid
    x_recon_bad = reconstruct_signal(t_bad, x_bad, fs_bad, t_recon)

    # Frequency domain: compute spectra of good and bad sampling
    # Use the sampled signals, zero-padded for resolution
    nfft = 2048
    freqs_good = np.fft.rfftfreq(nfft, 1 / fs_good)
    mag_good = 2.0 / len(x_good) * np.abs(np.fft.rfft(x_good, n=nfft))
    freqs_bad = np.fft.rfftfreq(nfft, 1 / fs_bad)
    mag_bad = 2.0 / len(x_bad) * np.abs(np.fft.rfft(x_bad, n=nfft))

    # --- Plot ---
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("Sampling & Aliasing", fontsize=14)

    # Panel 1: Original continuous signal
    ax = axes[0, 0]
    ax.plot(t_cont, x_cont, color="tab:blue", linewidth=1.5)
    ax.set_title(f"Original Continuous Signal ({freq} Hz)")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Amplitude")
    ax.grid(True, alpha=0.3)

    # Panel 2: Properly sampled (fs > 2*f)
    ax = axes[0, 1]
    ax.plot(t_cont, x_cont, color="tab:blue", alpha=0.4, linewidth=1,
            label="Continuous")
    ax.stem(t_good, x_good, linefmt="tab:green", markerfmt="go",
            basefmt="none", label=f"Samples (fs={fs_good} Hz)")
    ax.set_title(f"Properly Sampled (fs={fs_good} Hz > 2×{freq} Hz)")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Amplitude")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel 3: Undersampled — aliasing visible
    ax = axes[1, 0]
    ax.plot(t_cont, x_cont, color="tab:blue", alpha=0.4, linewidth=1,
            label=f"Original ({freq} Hz)")
    ax.stem(t_bad, x_bad, linefmt="tab:red", markerfmt="ro",
            basefmt="none", label=f"Samples (fs={fs_bad} Hz)")
    ax.plot(t_recon, x_recon_bad, color="tab:red", linewidth=1.5,
            linestyle="--", alpha=0.8,
            label=f"Reconstructed ({aliased_freq} Hz alias)")
    ax.set_title(f"Undersampled (fs={fs_bad} Hz < 2×{freq} Hz) — Aliasing!")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Amplitude")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel 4: Frequency domain comparison
    ax = axes[1, 1]
    ax.stem(freqs_good, mag_good, linefmt="tab:green", markerfmt="go",
            basefmt="none", label=f"Good sampling (fs={fs_good} Hz)")
    ax.stem(freqs_bad, mag_bad, linefmt="tab:red", markerfmt="r^",
            basefmt="none", label=f"Undersampled (fs={fs_bad} Hz)")
    ax.axvline(freq, color="tab:blue", linestyle=":", alpha=0.6,
               label=f"True freq ({freq} Hz)")
    if aliased_freq != freq:
        ax.axvline(aliased_freq, color="tab:red", linestyle=":", alpha=0.6,
                   label=f"Aliased freq ({aliased_freq} Hz)")
    ax.set_title("Frequency Domain — Where Energy Appears")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Magnitude")
    ax.set_xlim(0, max(fs_good, fs_bad) / 2 + 2)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Interactive widget (for Jupyter notebook)
# ---------------------------------------------------------------------------

def interactive_sampling():
    """Launch ipywidgets sliders for exploring sampling and aliasing."""
    try:
        from ipywidgets import interact, FloatSlider
    except ImportError:
        print("ipywidgets not available — run in Jupyter for interactive mode")
        return

    def _update(freq=5.0, fs=40.0):
        duration = 1.0
        t_cont, x_cont, t_sampled, x_sampled = sample_signal(
            freq, duration, fs
        )
        aliased_freq, is_aliased = demonstrate_aliasing(freq, fs)

        # Reconstruct from samples
        x_recon = reconstruct_signal(t_sampled, x_sampled, fs, t_cont)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

        # Time domain
        ax1.plot(t_cont, x_cont, color="tab:blue", alpha=0.4, linewidth=1,
                 label=f"Original ({freq} Hz)")
        ax1.stem(t_sampled, x_sampled, linefmt="tab:green", markerfmt="go",
                 basefmt="none", label=f"Samples (fs={fs} Hz)")
        ax1.plot(t_cont, x_recon, color="tab:red", linewidth=1.5,
                 linestyle="--", alpha=0.7, label="Reconstructed")

        status = "ALIASED" if is_aliased else "OK"
        color = "red" if is_aliased else "green"
        ax1.set_title(f"Sampling: {status} "
                      f"(Nyquist = {fs/2} Hz, fmax = {freq} Hz)",
                      color=color)
        ax1.set_xlabel("Time (s)")
        ax1.set_ylabel("Amplitude")
        ax1.legend(fontsize=8)
        ax1.grid(True, alpha=0.3)

        # Frequency domain
        nfft = 2048
        freqs = np.fft.rfftfreq(nfft, 1 / fs)
        mag = 2.0 / len(x_sampled) * np.abs(np.fft.rfft(x_sampled, n=nfft))
        ax2.stem(freqs, mag, linefmt="tab:blue", markerfmt="o",
                 basefmt="none")
        ax2.axvline(freq, color="tab:blue", linestyle=":", alpha=0.5,
                    label=f"True freq ({freq} Hz)")
        if is_aliased:
            ax2.axvline(aliased_freq, color="tab:red", linestyle=":",
                        alpha=0.7,
                        label=f"Aliased to {aliased_freq:.1f} Hz")
        ax2.set_title("Frequency Domain")
        ax2.set_xlabel("Frequency (Hz)")
        ax2.set_ylabel("Magnitude")
        ax2.set_xlim(0, fs / 2 + 1)
        ax2.legend(fontsize=8)
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

    interact(
        _update,
        freq=FloatSlider(value=5.0, min=1.0, max=50.0, step=0.5,
                         description="Signal (Hz)"),
        fs=FloatSlider(value=40.0, min=2.0, max=100.0, step=1.0,
                       description="Sample Rate"),
    )


# ---------------------------------------------------------------------------
# Standalone execution
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    matplotlib.use("Agg")
    fig = plot_sampling_demo()
    fig.savefig("sampling_aliasing.png", dpi=150)
    print("Saved sampling_aliasing.png")

    plt.show()
