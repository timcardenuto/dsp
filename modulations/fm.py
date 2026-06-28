"""
Frequency Modulation (FM)

In AM (Amplitude Modulation), you encode information by changing how loud a
carrier wave is.  FM takes a completely different approach: the carrier stays
at constant loudness, but its *frequency* speeds up and slows down to follow
the message signal.

Why does this matter?  Most real-world noise — static, lightning, electrical
interference — shows up as random changes in *amplitude*.  Because FM carries
its information in frequency, not amplitude, it naturally ignores most of that
noise.  This is why FM radio sounds cleaner than AM radio.

Key concepts:

  Frequency deviation (delta_f)
      How far the carrier's frequency swings above and below its resting
      (center) frequency.  A message value of +1 pushes the frequency up by
      delta_f; a value of -1 pulls it down by delta_f.

  Modulation index (beta)
      beta = delta_f / f_message.  This single number captures "how dramatic"
      the FM is.  Small beta (<1) is narrowband FM; large beta (>>1) is
      wideband FM.  Commercial FM radio uses beta ~ 5.

  Carson's rule (bandwidth estimate)
      BW ≈ 2 * (delta_f + f_message).  This tells you how much spectrum an FM
      signal occupies.  Higher deviation or higher message frequency → wider
      bandwidth.

  Instantaneous frequency
      At any instant the FM signal's frequency equals:
          f_inst(t) = f_carrier + delta_f * message(t)
      So plotting f_inst over time should look just like the original message —
      that's how you recover it (demodulation).

The math:
  An FM signal is   x(t) = A * cos(2*pi*f_c*t + 2*pi*delta_f * integral(message) dt)
  The key insight is that *phase is the integral of frequency*, so to make the
  instantaneous frequency follow the message, we integrate the message and put
  it in the phase argument.
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.signal import hilbert, spectrogram as _spectrogram


# ---------------------------------------------------------------------------
# Core functions (importable by notebook)
# ---------------------------------------------------------------------------

def fm_modulate(message, carrier_freq, fs, freq_deviation, amplitude=1.0):
    """Generate an FM-modulated signal.

    Parameters
    ----------
    message : ndarray – normalized message signal (values in roughly [-1, 1])
    carrier_freq : float – carrier frequency in Hz
    fs : float – sample rate in samples/second
    freq_deviation : float – maximum frequency deviation in Hz (delta_f)
    amplitude : float – carrier amplitude

    Returns
    -------
    fm_signal : ndarray – FM-modulated waveform (same length as message)

    Notes
    -----
    The instantaneous frequency is  f_carrier + freq_deviation * message(t).
    Phase is the integral of frequency, so we cumsum the message and scale
    appropriately.
    """
    N = len(message)
    t = np.arange(N) / fs

    # Phase contributed by the message: integral of (2*pi*freq_deviation*message)
    # Using cumulative sum as a discrete approximation to the integral.
    message_integral = np.cumsum(message) / fs
    phase = 2 * np.pi * carrier_freq * t + 2 * np.pi * freq_deviation * message_integral

    fm_signal = amplitude * np.cos(phase)
    return fm_signal


def fm_demodulate(signal, fs):
    """Demodulate an FM signal by extracting its instantaneous frequency.

    This uses the analytic signal (via the Hilbert transform) to compute the
    instantaneous phase, then differentiates to get instantaneous frequency.
    The result is proportional to the original message.

    Parameters
    ----------
    signal : ndarray – FM-modulated signal
    fs : float – sample rate in samples/second

    Returns
    -------
    demodulated : ndarray – recovered message (proportional, not calibrated)
    """
    inst_freq = instantaneous_frequency(signal, fs)
    # Remove the DC component (carrier frequency) to recover the message shape
    demodulated = inst_freq - np.mean(inst_freq)
    return demodulated


def instantaneous_frequency(signal, fs):
    """Compute the instantaneous frequency of a signal.

    Uses the Hilbert transform to form the analytic signal, then computes
    frequency as the derivative of the unwrapped instantaneous phase.

    Parameters
    ----------
    signal : ndarray – input signal
    fs : float – sample rate in samples/second

    Returns
    -------
    inst_freq : ndarray – instantaneous frequency in Hz (same length as signal)
    """
    analytic = hilbert(signal)
    inst_phase = np.unwrap(np.angle(analytic))
    # Frequency = d(phase)/dt / (2*pi).  Use gradient for same-length output.
    inst_freq = np.gradient(inst_phase, 1 / fs) / (2 * np.pi)
    return inst_freq


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def plot_fm_demo(carrier_freq=1000, message_freq=100, freq_deviation=200,
                 fs=20000, duration=0.05):
    """Generate a five-panel overview of FM modulation.

    Panels:
      1. Message signal (the information we want to transmit)
      2. Carrier signal (unmodulated)
      3. FM-modulated signal
      4. Instantaneous frequency of the FM signal (should track the message)
      5. Spectrogram showing frequency variation over time
    """
    t = np.arange(0, duration, 1 / fs)
    message = np.sin(2 * np.pi * message_freq * t)
    carrier = np.cos(2 * np.pi * carrier_freq * t)
    fm_signal = fm_modulate(message, carrier_freq, fs, freq_deviation)
    inst_freq = instantaneous_frequency(fm_signal, fs)

    beta = freq_deviation / message_freq
    bw = 2 * (freq_deviation + message_freq)

    fig, axes = plt.subplots(5, 1, figsize=(12, 14))
    fig.suptitle(
        f"Frequency Modulation (FM)\n"
        f"β = Δf / f_msg = {freq_deviation}/{message_freq} = {beta:.1f}    "
        f"Carson BW ≈ 2·(Δf + f_msg) = {bw:.0f} Hz",
        fontsize=13,
    )

    t_ms = t * 1000  # plot in milliseconds

    # Panel 1 — Message
    ax = axes[0]
    ax.plot(t_ms, message, color="tab:green")
    ax.set_title(f"Message Signal ({message_freq} Hz sine)")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude")
    ax.grid(True, alpha=0.3)

    # Panel 2 — Carrier
    ax = axes[1]
    ax.plot(t_ms, carrier, color="tab:blue", alpha=0.7)
    ax.set_title(f"Carrier ({carrier_freq} Hz, unmodulated)")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude")
    ax.grid(True, alpha=0.3)

    # Panel 3 — FM signal
    ax = axes[2]
    ax.plot(t_ms, fm_signal, color="tab:purple")
    ax.set_title("FM Modulated Signal (notice bunching / spreading of cycles)")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude")
    ax.grid(True, alpha=0.3)

    # Panel 4 — Instantaneous frequency
    ax = axes[3]
    ax.plot(t_ms, inst_freq, color="tab:red", label="Instantaneous freq")
    ax.axhline(carrier_freq, color="gray", linestyle="--", alpha=0.5,
               label=f"Carrier ({carrier_freq} Hz)")
    ax.axhline(carrier_freq + freq_deviation, color="gray", linestyle=":",
               alpha=0.4, label=f"± deviation ({freq_deviation} Hz)")
    ax.axhline(carrier_freq - freq_deviation, color="gray", linestyle=":",
               alpha=0.4)
    ax.set_title("Instantaneous Frequency (should track the message shape)")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Frequency (Hz)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel 5 — Spectrogram
    ax = axes[4]
    nperseg = min(256, len(fm_signal) // 4)
    f_spec, t_spec, Sxx = _spectrogram(fm_signal, fs, nperseg=nperseg)
    Sxx_db = 10 * np.log10(Sxx + 1e-12)
    im = ax.pcolormesh(t_spec * 1000, f_spec, Sxx_db, shading="gouraud",
                       cmap="viridis")
    ax.set_ylim(0, carrier_freq + 3 * freq_deviation)
    ax.set_title("Spectrogram (frequency content over time)")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Frequency (Hz)")
    fig.colorbar(im, ax=ax, label="Power (dB)")

    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Interactive widget (for Jupyter notebook)
# ---------------------------------------------------------------------------

def interactive_fm():
    """Launch ipywidgets sliders for exploring FM modulation.

    Sliders control carrier frequency, message frequency, and frequency
    deviation.  The plot updates live so you can see how each parameter
    changes the modulated signal and its instantaneous frequency.
    """
    try:
        from ipywidgets import interact, FloatSlider
    except ImportError:
        print("ipywidgets not available — run in Jupyter for interactive mode")
        return

    def _update(carrier_freq=1000, message_freq=100, freq_deviation=200):
        fs = 20000
        duration = 0.05
        t = np.arange(0, duration, 1 / fs)
        message = np.sin(2 * np.pi * message_freq * t)
        fm_signal = fm_modulate(message, carrier_freq, fs, freq_deviation)
        inst_freq = instantaneous_frequency(fm_signal, fs)

        beta = freq_deviation / message_freq if message_freq > 0 else 0
        bw = 2 * (freq_deviation + message_freq)

        fig, axes = plt.subplots(3, 1, figsize=(12, 8))
        fig.suptitle(
            f"FM Interactive  |  β = {beta:.2f}  |  Carson BW ≈ {bw:.0f} Hz",
            fontsize=12,
        )
        t_ms = t * 1000

        axes[0].plot(t_ms, message, color="tab:green")
        axes[0].set_title("Message")
        axes[0].set_ylabel("Amplitude")
        axes[0].grid(True, alpha=0.3)

        axes[1].plot(t_ms, fm_signal, color="tab:purple")
        axes[1].set_title("FM Signal")
        axes[1].set_ylabel("Amplitude")
        axes[1].grid(True, alpha=0.3)

        axes[2].plot(t_ms, inst_freq, color="tab:red")
        axes[2].axhline(carrier_freq, color="gray", linestyle="--", alpha=0.5)
        axes[2].set_title("Instantaneous Frequency")
        axes[2].set_xlabel("Time (ms)")
        axes[2].set_ylabel("Freq (Hz)")
        axes[2].grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

    interact(
        _update,
        carrier_freq=FloatSlider(
            value=1000, min=200, max=5000, step=100,
            description="Carrier (Hz)",
        ),
        message_freq=FloatSlider(
            value=100, min=10, max=500, step=10,
            description="Message (Hz)",
        ),
        freq_deviation=FloatSlider(
            value=200, min=10, max=2000, step=10,
            description="Deviation (Hz)",
        ),
    )


# ---------------------------------------------------------------------------
# Standalone execution
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    matplotlib.use("Agg")

    fig = plot_fm_demo()
    fig.savefig("fm_modulation.png", dpi=150)
    print("Saved fm_modulation.png")

    plt.show()
