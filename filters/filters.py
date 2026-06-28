"""
Filters — Selecting Which Frequencies to Keep

A filter is like a gatekeeper for frequencies. Every real-world signal is a
mixture of many frequencies (as we saw with the FFT). A filter lets you choose
which frequencies pass through and which get blocked. This is essential because
signals almost always contain unwanted content — noise, interference from other
sources, or frequencies outside the range you care about.

There are three basic filter types:
  - Low-pass: keeps low frequencies, blocks high ones.
    Think of it as "keep the slow changes, remove the fast wiggles."
  - High-pass: keeps high frequencies, blocks low ones.
    Think of it as "remove the slow drift, keep the fast details."
  - Band-pass: keeps a specific range of frequencies, blocks everything else.
    Think of it as "tune in to one station and ignore the rest."

This module uses FIR (Finite Impulse Response) filters. "Finite" means the
filter only looks at a fixed number of past samples to compute each output
sample. The number of past samples is called the filter "order" — higher order
gives a sharper cutoff between passed and blocked frequencies, but costs more
computation.

The MATLAB files bandpass_filter.m and ResolutionBandwidth.m demonstrated
filtering by averaging samples at specific intervals — a manual approach to
frequency selection. Here we use scipy.signal.firwin, which designs the filter
coefficients mathematically for precise control over which frequencies pass.
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.signal import firwin, lfilter, freqz


# ---------------------------------------------------------------------------
# Core functions (importable by notebook)
# ---------------------------------------------------------------------------

def design_lowpass(cutoff, fs, order=51):
    """Design a low-pass FIR filter.

    A low-pass filter passes frequencies below `cutoff` and attenuates
    frequencies above it.

    Parameters
    ----------
    cutoff : float – cutoff frequency in Hz
    fs : float – sample rate in samples/second
    order : int – filter order (number of taps); higher = sharper cutoff

    Returns
    -------
    coeffs : ndarray – FIR filter coefficients
    """
    # firwin expects the cutoff as a fraction of the Nyquist frequency,
    # but passing fs= lets us specify it directly in Hz.
    coeffs = firwin(numtaps=order, cutoff=cutoff, fs=fs, pass_zero="lowpass")
    return coeffs


def design_highpass(cutoff, fs, order=51):
    """Design a high-pass FIR filter.

    A high-pass filter passes frequencies above `cutoff` and attenuates
    frequencies below it.

    Parameters
    ----------
    cutoff : float – cutoff frequency in Hz
    fs : float – sample rate in samples/second
    order : int – filter order (must be odd for high-pass FIR)

    Returns
    -------
    coeffs : ndarray – FIR filter coefficients
    """
    # High-pass FIR filters require an odd number of taps.
    if order % 2 == 0:
        order += 1
    coeffs = firwin(numtaps=order, cutoff=cutoff, fs=fs, pass_zero="highpass")
    return coeffs


def design_bandpass(low_cutoff, high_cutoff, fs, order=51):
    """Design a band-pass FIR filter.

    A band-pass filter passes frequencies between `low_cutoff` and
    `high_cutoff` and attenuates everything outside that range.

    Parameters
    ----------
    low_cutoff : float – lower cutoff frequency in Hz
    high_cutoff : float – upper cutoff frequency in Hz
    fs : float – sample rate in samples/second
    order : int – filter order (must be odd for band-pass FIR)

    Returns
    -------
    coeffs : ndarray – FIR filter coefficients
    """
    if order % 2 == 0:
        order += 1
    coeffs = firwin(
        numtaps=order,
        cutoff=[low_cutoff, high_cutoff],
        fs=fs,
        pass_zero="bandpass",
    )
    return coeffs


def apply_filter(signal, coeffs):
    """Apply a FIR filter to a signal.

    This convolves the filter coefficients with the signal. Each output
    sample is a weighted sum of nearby input samples, where the weights
    are the filter coefficients.

    Parameters
    ----------
    signal : ndarray – input signal samples
    coeffs : ndarray – FIR filter coefficients (from design_* functions)

    Returns
    -------
    filtered : ndarray – filtered signal (same length as input)
    """
    # lfilter applies the filter in the time domain.  For an FIR filter the
    # denominator polynomial is just [1] (no feedback).
    filtered = lfilter(coeffs, 1.0, signal)
    return filtered


def plot_frequency_response(coeffs, fs, title="Filter Frequency Response",
                            ax=None):
    """Plot the magnitude response of a filter.

    The frequency response shows how much the filter amplifies or attenuates
    each frequency.  A value of 0 dB means the frequency passes unchanged;
    -20 dB means it is reduced to 1/10th of its original amplitude.

    Parameters
    ----------
    coeffs : ndarray – FIR filter coefficients
    fs : float – sample rate in samples/second
    title : str – plot title
    ax : matplotlib Axes or None – axes to plot on (creates new figure if None)

    Returns
    -------
    fig : Figure or None – the figure (None if ax was provided)
    """
    w, h = freqz(coeffs, worN=2048, fs=fs)
    magnitude_db = 20 * np.log10(np.abs(h) + 1e-12)

    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 4))

    ax.plot(w, magnitude_db, color="tab:blue")
    ax.set_title(title)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Magnitude (dB)")
    ax.set_ylim(-80, 5)
    ax.grid(True, alpha=0.3)
    ax.axhline(-3, color="tab:red", linestyle="--", alpha=0.5, label="-3 dB")
    ax.legend(fontsize=8)

    if fig is not None:
        plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def _compute_spectrum(signal, fs):
    """Compute single-sided amplitude spectrum (internal helper)."""
    N = len(signal)
    Y = np.fft.fft(signal)
    freqs = np.fft.fftfreq(N, 1 / fs)[:N // 2]
    magnitude = 2.0 / N * np.abs(Y[:N // 2])
    return freqs, magnitude


def plot_filter_demo(fs=8000, duration=0.1):
    """Generate a four-panel figure showing what filtering does.

    Scenario: we have a 300 Hz signal we care about, contaminated by a
    1500 Hz tone and random noise.  A low-pass filter at 600 Hz removes
    the unwanted high-frequency content while preserving the desired signal.
    """
    # -- Build a noisy signal ------------------------------------------------
    desired_freq = 300      # Hz — the signal we want to keep
    unwanted_freq = 1500    # Hz — interference we want to remove
    t = np.arange(0, duration, 1 / fs)

    desired = 1.0 * np.sin(2 * np.pi * desired_freq * t)
    unwanted = 0.5 * np.sin(2 * np.pi * unwanted_freq * t)
    noise = 0.2 * np.random.default_rng(42).standard_normal(len(t))
    noisy_signal = desired + unwanted + noise

    # -- Design and apply a low-pass filter ----------------------------------
    cutoff = 600  # Hz
    order = 101
    coeffs = design_lowpass(cutoff, fs, order=order)
    filtered_signal = apply_filter(noisy_signal, coeffs)

    # -- Compute spectra -----------------------------------------------------
    freqs_orig, mag_orig = _compute_spectrum(noisy_signal, fs)
    freqs_filt, mag_filt = _compute_spectrum(filtered_signal, fs)

    # -- Frequency response of the filter ------------------------------------
    w, h = freqz(coeffs, worN=2048, fs=fs)
    mag_resp_db = 20 * np.log10(np.abs(h) + 1e-12)

    # -- Plot ----------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle(
        "Filtering Demo — Low-pass at 600 Hz removes 1500 Hz interference",
        fontsize=13,
    )

    # Panel 1: original noisy signal in time domain
    ax = axes[0, 0]
    ax.plot(t * 1000, noisy_signal, color="tab:blue", alpha=0.8)
    ax.set_title("Original Signal (300 Hz + 1500 Hz + noise)")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude")
    ax.grid(True, alpha=0.3)

    # Panel 2: filter frequency response
    ax = axes[0, 1]
    ax.plot(w, mag_resp_db, color="tab:orange")
    ax.axvline(cutoff, color="tab:red", linestyle="--", alpha=0.6,
               label=f"Cutoff = {cutoff} Hz")
    ax.axvline(desired_freq, color="tab:green", linestyle=":", alpha=0.6,
               label=f"Desired = {desired_freq} Hz")
    ax.axvline(unwanted_freq, color="tab:purple", linestyle=":", alpha=0.6,
               label=f"Unwanted = {unwanted_freq} Hz")
    ax.set_title("Low-pass Filter Frequency Response")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Magnitude (dB)")
    ax.set_ylim(-80, 5)
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    # Panel 3: filtered signal in time domain
    ax = axes[1, 0]
    ax.plot(t * 1000, desired, color="tab:green", alpha=0.4,
            label="Pure 300 Hz (reference)")
    ax.plot(t * 1000, filtered_signal, color="tab:red", alpha=0.8,
            label="Filtered output")
    ax.set_title("Filtered Signal vs Pure 300 Hz Reference")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel 4: before/after frequency spectra
    ax = axes[1, 1]
    ax.plot(freqs_orig, mag_orig, color="tab:blue", alpha=0.6,
            label="Before filtering")
    ax.plot(freqs_filt, mag_filt, color="tab:red", alpha=0.8,
            label="After filtering")
    ax.set_title("Frequency Spectra — Before vs After")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Amplitude")
    ax.set_xlim(0, 2500)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Interactive widget (for Jupyter notebook)
# ---------------------------------------------------------------------------

def interactive_filter():
    """Launch ipywidgets sliders for exploring filter behavior.

    Sliders let you adjust:
      - cutoff frequency: where the filter starts blocking
      - filter order: how sharp the transition is
      - signal frequencies: choose which tones are in the input
    """
    try:
        from ipywidgets import interact, IntSlider, FloatSlider
    except ImportError:
        print("ipywidgets not available — run in Jupyter for interactive mode")
        return

    def _update(cutoff=600, order=51, f_desired=300, f_unwanted=1500,
                noise_level=0.2):
        fs = 8000
        duration = 0.1
        t = np.arange(0, duration, 1 / fs)

        # Ensure order is odd (required for high-pass / band-pass FIR)
        if order % 2 == 0:
            order += 1

        # Build signal
        desired = np.sin(2 * np.pi * f_desired * t)
        unwanted = 0.5 * np.sin(2 * np.pi * f_unwanted * t)
        noise = noise_level * np.random.default_rng(42).standard_normal(len(t))
        signal = desired + unwanted + noise

        # Design and apply filter
        coeffs = design_lowpass(cutoff, fs, order=order)
        filtered = apply_filter(signal, coeffs)

        # Spectra
        freqs_orig, mag_orig = _compute_spectrum(signal, fs)
        freqs_filt, mag_filt = _compute_spectrum(filtered, fs)

        # Frequency response
        w, h = freqz(coeffs, worN=2048, fs=fs)
        mag_db = 20 * np.log10(np.abs(h) + 1e-12)

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        fig.suptitle("Interactive Filter Explorer", fontsize=13)

        # Time domain — original
        ax = axes[0, 0]
        ax.plot(t * 1000, signal, color="tab:blue", alpha=0.8)
        ax.set_title(f"Original ({f_desired} Hz + {f_unwanted} Hz + noise)")
        ax.set_xlabel("Time (ms)")
        ax.set_ylabel("Amplitude")
        ax.grid(True, alpha=0.3)

        # Filter response
        ax = axes[0, 1]
        ax.plot(w, mag_db, color="tab:orange")
        ax.axvline(cutoff, color="tab:red", linestyle="--", alpha=0.6,
                   label=f"Cutoff = {cutoff} Hz")
        ax.set_title(f"Low-pass Response (order={order})")
        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel("Magnitude (dB)")
        ax.set_ylim(-80, 5)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        # Time domain — filtered
        ax = axes[1, 0]
        ax.plot(t * 1000, filtered, color="tab:red", alpha=0.8)
        ax.set_title("Filtered Signal")
        ax.set_xlabel("Time (ms)")
        ax.set_ylabel("Amplitude")
        ax.grid(True, alpha=0.3)

        # Spectra comparison
        ax = axes[1, 1]
        ax.plot(freqs_orig, mag_orig, color="tab:blue", alpha=0.6,
                label="Before")
        ax.plot(freqs_filt, mag_filt, color="tab:red", alpha=0.8,
                label="After")
        ax.set_title("Frequency Spectra — Before vs After")
        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel("Amplitude")
        ax.set_xlim(0, fs / 2)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

    interact(
        _update,
        cutoff=FloatSlider(
            value=600, min=50, max=3500, step=50,
            description="Cutoff (Hz)",
        ),
        order=IntSlider(
            value=51, min=11, max=201, step=10,
            description="Filter Order",
        ),
        f_desired=FloatSlider(
            value=300, min=50, max=2000, step=50,
            description="Desired (Hz)",
        ),
        f_unwanted=FloatSlider(
            value=1500, min=50, max=3500, step=50,
            description="Unwanted (Hz)",
        ),
        noise_level=FloatSlider(
            value=0.2, min=0.0, max=1.0, step=0.05,
            description="Noise Level",
        ),
    )


# ---------------------------------------------------------------------------
# Standalone execution
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    matplotlib.use("Agg")

    fig1 = plot_filter_demo()
    fig1.savefig("filter_demo.png", dpi=150)
    print("Saved filter_demo.png")

    # Also save individual frequency responses for the three filter types
    fs = 8000
    fig2, axes = plt.subplots(1, 3, figsize=(15, 4))
    fig2.suptitle("Filter Types — Frequency Responses", fontsize=13)

    lp = design_lowpass(1000, fs, order=51)
    w, h = freqz(lp, worN=2048, fs=fs)
    axes[0].plot(w, 20 * np.log10(np.abs(h) + 1e-12), color="tab:blue")
    axes[0].set_title("Low-pass (cutoff = 1000 Hz)")
    axes[0].set_xlabel("Frequency (Hz)")
    axes[0].set_ylabel("Magnitude (dB)")
    axes[0].set_ylim(-80, 5)
    axes[0].grid(True, alpha=0.3)

    hp = design_highpass(1000, fs, order=51)
    w, h = freqz(hp, worN=2048, fs=fs)
    axes[1].plot(w, 20 * np.log10(np.abs(h) + 1e-12), color="tab:orange")
    axes[1].set_title("High-pass (cutoff = 1000 Hz)")
    axes[1].set_xlabel("Frequency (Hz)")
    axes[1].set_ylabel("Magnitude (dB)")
    axes[1].set_ylim(-80, 5)
    axes[1].grid(True, alpha=0.3)

    bp = design_bandpass(800, 1200, fs, order=101)
    w, h = freqz(bp, worN=2048, fs=fs)
    axes[2].plot(w, 20 * np.log10(np.abs(h) + 1e-12), color="tab:green")
    axes[2].set_title("Band-pass (800-1200 Hz)")
    axes[2].set_xlabel("Frequency (Hz)")
    axes[2].set_ylabel("Magnitude (dB)")
    axes[2].set_ylim(-80, 5)
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    fig2.savefig("filter_types.png", dpi=150)
    print("Saved filter_types.png")

    plt.show()
