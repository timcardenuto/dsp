"""
Amplitude Modulation (AM)

Imagine you want to send music from a radio station to your car stereo. The
audio signal (music, voice) has low frequencies — roughly 20 Hz to 20 kHz —
that can't travel far through the air on their own. AM radio solves this by
"piggy-backing" the audio onto a high-frequency wave that *can* travel long
distances. That high-frequency wave is called the **carrier**, and the audio
is called the **message** (or baseband) signal.

How it works:
  1. Start with a carrier wave — a pure sine wave at a high frequency (e.g.
     540 kHz for AM radio). By itself it carries no information.
  2. Vary (modulate) the carrier's *amplitude* in proportion to the message
     signal. When the message is loud, the carrier gets taller; when the
     message is quiet, the carrier gets shorter.
  3. The result is a wave whose outline (envelope) traces out the original
     message. A receiver can recover the audio by extracting that envelope.

Key terms:
  - Carrier wave: the high-frequency sine wave used to transport the signal.
  - Message signal: the information you actually want to send.
  - Modulation index (m): controls how much the carrier amplitude changes.
        m = (peak message amplitude) / (carrier amplitude).
    At m = 1.0 the carrier amplitude swings between 0 and 2x its original
    value. At m < 1.0 the envelope never hits zero. At m > 1.0 the signal
    is "overmodulated" — the envelope crosses zero and distorts, making
    simple demodulation fail.
  - Envelope detection: the simplest demodulation method. Rectify the signal
    (take the absolute value — like a diode in a crystal radio) then smooth
    it with a low-pass filter to strip out the carrier and recover the
    message. This is exactly how the cheapest AM radios work.

The math:
  carrier(t)   = Ac * sin(2*pi*fc*t)
  message(t)   = Am * sin(2*pi*fm*t)
  AM signal(t) = Ac * [1 + m * message(t)/Am] * sin(2*pi*fc*t)

In the frequency domain, AM produces three components:
  - The carrier at fc
  - An upper sideband at fc + fm
  - A lower sideband at fc - fm
This is why AM stations need a channel bandwidth of 2x the highest audio
frequency they transmit.
"""

import numpy as np
from scipy.signal import butter, filtfilt
import matplotlib
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Core functions (importable by notebook)
# ---------------------------------------------------------------------------

def am_modulate(message, carrier_freq, fs, mod_index=1.0, carrier_amp=1.0):
    """Create an AM-modulated signal.

    The modulation formula is:
        s(t) = Ac * [1 + m * message_normalized(t)] * sin(2*pi*fc*t)

    where message_normalized is the message scaled to the range [-1, 1].

    Parameters
    ----------
    message : ndarray – message (baseband) signal samples
    carrier_freq : float – carrier frequency in Hz
    fs : float – sample rate in samples/second
    mod_index : float – modulation index (0 to 1 for normal AM; >1 overmodulates)
    carrier_amp : float – peak amplitude of the unmodulated carrier

    Returns
    -------
    t : ndarray – time vector
    modulated : ndarray – AM-modulated signal
    carrier : ndarray – the bare carrier (useful for plotting)
    """
    N = len(message)
    t = np.arange(N) / fs
    carrier = carrier_amp * np.sin(2 * np.pi * carrier_freq * t)

    # Normalize message to [-1, 1] so mod_index directly controls depth
    msg_peak = np.max(np.abs(message))
    if msg_peak > 0:
        message_norm = message / msg_peak
    else:
        message_norm = message

    modulated = carrier_amp * (1 + mod_index * message_norm) * np.sin(2 * np.pi * carrier_freq * t)
    return t, modulated, carrier


def am_demodulate(signal, fs, carrier_freq):
    """Recover the message from an AM signal via envelope detection.

    This mimics how a simple crystal/diode radio works:
      1. Rectify — take the absolute value (like a diode that blocks the
         negative half of each cycle).
      2. Low-pass filter — smooth out the rapid carrier oscillations to
         reveal the slow-moving envelope, which is the original message.

    The cutoff frequency is set midway between the message bandwidth and the
    carrier, so it passes the message but rejects the carrier.

    Parameters
    ----------
    signal : ndarray – AM-modulated signal
    fs : float – sample rate
    carrier_freq : float – carrier frequency (used to design the filter)

    Returns
    -------
    envelope : ndarray – recovered message (envelope of the AM signal)
    """
    # Step 1: Rectify (absolute value, like a diode)
    rectified = np.abs(signal)

    # Step 2: Low-pass filter to strip out the carrier frequency ripple.
    # Cutoff at carrier_freq / 2 keeps the message and rejects the carrier.
    cutoff = carrier_freq / 2
    nyquist = fs / 2
    if cutoff >= nyquist:
        # Fallback: simple moving average (matches DiodeRadioSim.m approach)
        window = max(1, int(fs / carrier_freq))
        kernel = np.ones(window) / window
        envelope = np.convolve(rectified, kernel, mode="same")
    else:
        b, a = butter(4, cutoff / nyquist, btype="low")
        envelope = filtfilt(b, a, rectified)

    # Remove DC offset — the envelope sits on top of a DC bias from the
    # carrier. Subtract the mean so the recovered signal is centered at zero,
    # matching the original message shape.
    envelope = envelope - np.mean(envelope)

    return envelope


def generate_am_signal(msg_freq=10, carrier_freq=100, fs=10000, duration=0.3,
                       mod_index=1.0, msg_amp=1.0):
    """Generate a complete AM example: message, carrier, and modulated signal.

    This is a convenience function that creates a simple sinusoidal message
    and modulates it onto a carrier.  Useful for quick experiments.

    Parameters
    ----------
    msg_freq : float – message frequency in Hz (default 10)
    carrier_freq : float – carrier frequency in Hz (default 100)
    fs : float – sample rate (default 10000)
    duration : float – signal duration in seconds (default 0.3)
    mod_index : float – modulation index (default 1.0)
    msg_amp : float – message amplitude (default 1.0)

    Returns
    -------
    t : ndarray – time vector
    message : ndarray – the message signal
    carrier : ndarray – the bare carrier
    modulated : ndarray – the AM signal
    """
    t = np.arange(0, duration, 1 / fs)
    message = msg_amp * np.sin(2 * np.pi * msg_freq * t)
    _, modulated, carrier = am_modulate(message, carrier_freq, fs, mod_index)
    return t, message, carrier, modulated


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_am_demo(msg_freq=10, carrier_freq=100, fs=10000, duration=0.3,
                 mod_index=0.8):
    """Generate a five-panel figure showing the full AM story.

    Panels:
      1. Message signal (the audio you want to send)
      2. Carrier signal (the high-frequency "taxi")
      3. AM modulated signal with envelope overlay
      4. Frequency spectrum showing carrier and sidebands
      5. Demodulated (recovered) signal vs original message
    """
    t, message, carrier, modulated = generate_am_signal(
        msg_freq=msg_freq, carrier_freq=carrier_freq, fs=fs,
        duration=duration, mod_index=mod_index,
    )

    # Demodulate
    recovered = am_demodulate(modulated, fs, carrier_freq)

    # Compute spectrum of the AM signal
    N = len(modulated)
    freqs = np.fft.fftfreq(N, 1 / fs)[:N // 2]
    spectrum = 2.0 / N * np.abs(np.fft.fft(modulated)[:N // 2])

    # Compute envelope for overlay on panel 3
    upper_env = np.abs(modulated)
    cutoff = carrier_freq / 2
    nyquist = fs / 2
    if cutoff < nyquist:
        b, a = butter(4, cutoff / nyquist, btype="low")
        upper_env = filtfilt(b, a, np.abs(modulated))

    fig, axes = plt.subplots(5, 1, figsize=(12, 14))
    fig.suptitle("Amplitude Modulation (AM) — From Message to Recovery",
                 fontsize=14, fontweight="bold")
    t_ms = t * 1000  # plot in milliseconds

    # Panel 1 — Message signal
    ax = axes[0]
    ax.plot(t_ms, message, color="tab:blue")
    ax.set_title(f"Message Signal ({msg_freq} Hz) — the information to send")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude")
    ax.grid(True, alpha=0.3)

    # Panel 2 — Carrier
    ax = axes[1]
    ax.plot(t_ms, carrier, color="tab:orange", linewidth=0.5)
    ax.set_title(f"Carrier Signal ({carrier_freq} Hz) — the high-frequency transport wave")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude")
    ax.grid(True, alpha=0.3)

    # Panel 3 — AM signal with envelope
    ax = axes[2]
    ax.plot(t_ms, modulated, color="tab:green", linewidth=0.5, label="AM signal")
    ax.plot(t_ms, upper_env, color="tab:red", linewidth=2, label="Envelope (message shape)")
    ax.plot(t_ms, -upper_env, color="tab:red", linewidth=2)
    ax.set_title(f"AM Modulated Signal (mod index = {mod_index})")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel 4 — Frequency spectrum
    ax = axes[3]
    ax.plot(freqs, spectrum, color="tab:purple")
    # Zoom in around the carrier
    view_width = carrier_freq * 2
    ax.set_xlim(0, view_width)
    ax.axvline(carrier_freq, color="tab:orange", linestyle="--", alpha=0.7,
               label=f"Carrier ({carrier_freq} Hz)")
    ax.axvline(carrier_freq - msg_freq, color="tab:blue", linestyle="--",
               alpha=0.7, label=f"Lower sideband ({carrier_freq - msg_freq} Hz)")
    ax.axvline(carrier_freq + msg_freq, color="tab:blue", linestyle="--",
               alpha=0.7, label=f"Upper sideband ({carrier_freq + msg_freq} Hz)")
    ax.set_title("Frequency Spectrum — carrier plus sidebands")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Amplitude")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel 5 — Demodulated signal
    ax = axes[4]
    # Scale recovered to match message amplitude for comparison
    msg_peak = np.max(np.abs(message))
    rec_peak = np.max(np.abs(recovered))
    if rec_peak > 0:
        recovered_scaled = recovered * (msg_peak / rec_peak)
    else:
        recovered_scaled = recovered
    ax.plot(t_ms, message, color="tab:blue", alpha=0.5, label="Original message")
    ax.plot(t_ms, recovered_scaled, color="tab:red", label="Recovered (envelope detection)")
    ax.set_title("Demodulated Signal — recovered by envelope detection")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Interactive widget (for Jupyter notebook)
# ---------------------------------------------------------------------------

def interactive_am():
    """Launch ipywidgets sliders for exploring AM modulation.

    Sliders:
      - Carrier frequency: see how the "taxi" wave speed changes
      - Message frequency: change the information signal
      - Modulation index: watch the envelope depth change; push past 1.0
        to see overmodulation distortion
    """
    try:
        from ipywidgets import interact, FloatSlider
    except ImportError:
        print("ipywidgets not available — run in Jupyter for interactive mode")
        return

    def _update(carrier_freq=100, msg_freq=10, mod_index=0.8):
        fs = 10000
        duration = 0.3
        t, message, carrier, modulated = generate_am_signal(
            msg_freq=msg_freq, carrier_freq=carrier_freq, fs=fs,
            duration=duration, mod_index=mod_index,
        )
        recovered = am_demodulate(modulated, fs, carrier_freq)

        # Scale recovered for visual comparison
        msg_peak = np.max(np.abs(message))
        rec_peak = np.max(np.abs(recovered))
        if rec_peak > 0:
            recovered_scaled = recovered * (msg_peak / rec_peak)
        else:
            recovered_scaled = recovered

        t_ms = t * 1000

        fig, axes = plt.subplots(3, 1, figsize=(12, 8))
        fig.suptitle(f"AM Explorer — mod index = {mod_index:.2f}"
                     + (" (OVERMODULATED!)" if mod_index > 1.0 else ""),
                     fontsize=13, fontweight="bold",
                     color="red" if mod_index > 1.0 else "black")

        # Compute envelope for overlay
        upper_env = np.abs(modulated)
        cutoff = carrier_freq / 2
        nyquist = fs / 2
        if cutoff < nyquist:
            b, a = butter(4, cutoff / nyquist, btype="low")
            upper_env = filtfilt(b, a, np.abs(modulated))

        axes[0].plot(t_ms, message, color="tab:blue")
        axes[0].set_title(f"Message ({msg_freq:.0f} Hz)")
        axes[0].set_ylabel("Amplitude")
        axes[0].grid(True, alpha=0.3)

        axes[1].plot(t_ms, modulated, color="tab:green", linewidth=0.5)
        axes[1].plot(t_ms, upper_env, color="tab:red", linewidth=2)
        axes[1].plot(t_ms, -upper_env, color="tab:red", linewidth=2)
        axes[1].set_title(f"AM Signal (carrier {carrier_freq:.0f} Hz)")
        axes[1].set_ylabel("Amplitude")
        axes[1].grid(True, alpha=0.3)

        axes[2].plot(t_ms, message, color="tab:blue", alpha=0.5, label="Original")
        axes[2].plot(t_ms, recovered_scaled, color="tab:red", label="Recovered")
        axes[2].set_title("Demodulated")
        axes[2].set_xlabel("Time (ms)")
        axes[2].set_ylabel("Amplitude")
        axes[2].legend(fontsize=8)
        axes[2].grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

    interact(
        _update,
        carrier_freq=FloatSlider(value=100, min=50, max=500, step=10,
                                 description="Carrier (Hz)"),
        msg_freq=FloatSlider(value=10, min=1, max=50, step=1,
                             description="Message (Hz)"),
        mod_index=FloatSlider(value=0.8, min=0.1, max=2.0, step=0.1,
                              description="Mod Index"),
    )


# ---------------------------------------------------------------------------
# Standalone execution
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    matplotlib.use("Agg")

    fig1 = plot_am_demo()
    fig1.savefig("am_modulation.png", dpi=150)
    print("Saved am_modulation.png")

    # Also generate an overmodulation comparison
    fig2 = plot_am_demo(mod_index=1.5)
    fig2.savefig("am_overmodulated.png", dpi=150)
    print("Saved am_overmodulated.png")

    plt.close("all")
