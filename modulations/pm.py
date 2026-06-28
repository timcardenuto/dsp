"""
Phase Modulation (PM)

Imagine a spinning wheel that completes one full rotation every cycle. At any
instant, the wheel's angle — its "phase" — tells you where it is in that cycle,
from 0 degrees (just starting) through 360 degrees (back to the start).

A plain carrier wave spins at a constant speed, so its phase increases at a
steady rate. Phase Modulation works by nudging that phase forward or backward
in proportion to a message signal. When the message is positive, the carrier's
phase jumps ahead; when the message is negative, it lags behind. The stronger
the message, the bigger the nudge.

Key parameters:
  - Carrier frequency (fc): the base spin rate of the unmodulated wave.
  - Phase deviation (kp): how many radians the phase shifts per unit of
    message amplitude. A larger kp means the phase swings more dramatically,
    which widens the bandwidth.

PM vs FM — close cousins:
  Phase Modulation and Frequency Modulation are mathematically siblings.
  Frequency is the *rate of change* of phase, so:
    - PM of a message m(t) produces the same waveform as FM of m'(t)
      (the derivative of the message).
    - FM of a message m(t) produces the same waveform as PM of the integral
      of m(t).
  In practice this means PM and FM look very similar, but PM tracks the
  message's *shape* in its phase, while FM tracks the message's shape in its
  *instantaneous frequency*.
"""

import numpy as np
from scipy.signal import hilbert
import matplotlib
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Core functions (importable by notebook)
# ---------------------------------------------------------------------------

def pm_modulate(message, carrier_freq, fs, phase_deviation):
    """Apply phase modulation to a message signal.

    The PM formula is:
        s(t) = cos(2*pi*fc*t + kp * m(t))

    where kp is the phase deviation constant and m(t) is the message.

    Parameters
    ----------
    message : ndarray – the baseband message signal (arbitrary units)
    carrier_freq : float – carrier frequency in Hz
    fs : float – sample rate in samples/second
    phase_deviation : float – kp, phase shift in radians per unit of message
                      amplitude. For example, kp = pi means the phase swings
                      +/- 180 degrees when the message hits +/- 1.

    Returns
    -------
    t : ndarray – time vector (seconds)
    modulated : ndarray – PM waveform
    """
    N = len(message)
    t = np.arange(N) / fs
    modulated = np.cos(2 * np.pi * carrier_freq * t + phase_deviation * message)
    return t, modulated


def pm_demodulate(signal, fs, carrier_freq):
    """Recover the message from a PM waveform using the analytic signal.

    Strategy: use the Hilbert transform to get the instantaneous phase of the
    received signal, then subtract the carrier's linearly-increasing phase to
    isolate the phase deviation that carries the message.

    Parameters
    ----------
    signal : ndarray – received PM waveform
    fs : float – sample rate in samples/second
    carrier_freq : float – carrier frequency in Hz (must match the modulator)

    Returns
    -------
    t : ndarray – time vector (seconds)
    recovered : ndarray – estimated message signal (scaled by kp)
    """
    N = len(signal)
    t = np.arange(N) / fs

    # Analytic signal via Hilbert transform
    analytic = hilbert(signal)
    instantaneous_phase = np.unwrap(np.angle(analytic))

    # The carrier contributes a linearly-growing phase: 2*pi*fc*t
    carrier_phase = 2 * np.pi * carrier_freq * t

    # What remains after subtracting the carrier is kp * m(t)
    recovered = instantaneous_phase - carrier_phase

    # Remove any DC offset introduced by the transform
    recovered = recovered - np.mean(recovered)

    return t, recovered


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def plot_pm_demo(carrier_freq=1000, message_freq=100, fs=50000,
                 phase_deviation=np.pi, duration=0.03):
    """Generate a five-panel figure illustrating Phase Modulation.

    Panels:
        1. Message signal
        2. Unmodulated carrier
        3. PM-modulated waveform
        4. Extracted phase vs time (shows the phase tracks the message)
        5. PM vs FM comparison (side-by-side overlay)
    """
    t = np.arange(0, duration, 1 / fs)
    message = np.cos(2 * np.pi * message_freq * t)

    # ---- Carrier ----
    carrier = np.cos(2 * np.pi * carrier_freq * t)

    # ---- PM ----
    _, pm_signal = pm_modulate(message, carrier_freq, fs, phase_deviation)

    # ---- Extract instantaneous phase from the PM signal ----
    analytic_pm = hilbert(pm_signal)
    inst_phase_pm = np.unwrap(np.angle(analytic_pm))
    carrier_phase = 2 * np.pi * carrier_freq * t
    phase_deviation_extracted = inst_phase_pm - carrier_phase
    phase_deviation_extracted -= np.mean(phase_deviation_extracted)

    # ---- FM signal for comparison ----
    # FM embeds the message in the *frequency*, not the phase.
    # s_FM(t) = cos(2*pi*fc*t + 2*pi*kf * integral(m(t)))
    freq_deviation = carrier_freq * 0.5  # kf chosen for visible effect
    integral_message = np.cumsum(message) / fs
    fm_signal = np.cos(2 * np.pi * carrier_freq * t
                       + 2 * np.pi * freq_deviation * integral_message)

    # ---- Figure ----
    t_ms = t * 1000  # plot in milliseconds

    fig, axes = plt.subplots(5, 1, figsize=(12, 14))
    fig.suptitle("Phase Modulation (PM) — Overview", fontsize=14)

    # Panel 1: Message
    ax = axes[0]
    ax.plot(t_ms, message, color="tab:green")
    ax.set_title("Message Signal  m(t) = cos(2*pi*Fm*t)")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude")
    ax.grid(True, alpha=0.3)

    # Panel 2: Carrier
    ax = axes[1]
    ax.plot(t_ms, carrier, color="tab:blue")
    ax.set_title(f"Carrier Signal  cos(2*pi*{carrier_freq}*t)")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude")
    ax.grid(True, alpha=0.3)

    # Panel 3: PM waveform
    ax = axes[2]
    ax.plot(t_ms, pm_signal, color="tab:orange")
    ax.set_title(f"PM Modulated Signal  (kp = {phase_deviation:.2f} rad)")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude")
    ax.grid(True, alpha=0.3)

    # Panel 4: Extracted phase
    ax = axes[3]
    ax.plot(t_ms, phase_deviation_extracted, color="tab:red",
            label="Extracted phase deviation")
    ax.plot(t_ms, phase_deviation * message, color="tab:green",
            linestyle="--", alpha=0.7, label="kp * m(t)  (expected)")
    ax.set_title("Phase vs Time — Phase Tracks the Message")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Phase (radians)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel 5: PM vs FM comparison
    ax = axes[4]
    ax.plot(t_ms, pm_signal, color="tab:orange", alpha=0.8, label="PM")
    ax.plot(t_ms, fm_signal, color="tab:purple", alpha=0.8, label="FM")
    ax.set_title("PM vs FM Comparison — Same Carrier & Message")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Interactive widget (for Jupyter notebook)
# ---------------------------------------------------------------------------

def interactive_pm():
    """Launch ipywidgets sliders to explore PM interactively.

    Sliders:
      - Carrier frequency (200 – 3000 Hz)
      - Message frequency (10 – 500 Hz)
      - Phase deviation kp (0 – 3*pi radians)
    """
    try:
        from ipywidgets import interact, FloatSlider
    except ImportError:
        print("ipywidgets not available — run in Jupyter for interactive mode")
        return

    def _update(carrier_freq=1000, message_freq=100, phase_deviation=3.14):
        fs = 50000
        duration = 0.03
        t = np.arange(0, duration, 1 / fs)
        message = np.cos(2 * np.pi * message_freq * t)

        _, pm_signal = pm_modulate(message, carrier_freq, fs, phase_deviation)
        carrier = np.cos(2 * np.pi * carrier_freq * t)

        t_ms = t * 1000

        fig, axes = plt.subplots(3, 1, figsize=(12, 8))

        axes[0].plot(t_ms, message, color="tab:green")
        axes[0].set_title("Message Signal")
        axes[0].set_xlabel("Time (ms)")
        axes[0].set_ylabel("Amplitude")
        axes[0].grid(True, alpha=0.3)

        axes[1].plot(t_ms, carrier, color="tab:blue")
        axes[1].set_title(f"Carrier  ({carrier_freq:.0f} Hz)")
        axes[1].set_xlabel("Time (ms)")
        axes[1].set_ylabel("Amplitude")
        axes[1].grid(True, alpha=0.3)

        axes[2].plot(t_ms, pm_signal, color="tab:orange")
        axes[2].set_title(f"PM Signal  (kp = {phase_deviation:.2f} rad)")
        axes[2].set_xlabel("Time (ms)")
        axes[2].set_ylabel("Amplitude")
        axes[2].grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

    interact(
        _update,
        carrier_freq=FloatSlider(
            value=1000, min=200, max=3000, step=50,
            description="Carrier (Hz)"),
        message_freq=FloatSlider(
            value=100, min=10, max=500, step=10,
            description="Message (Hz)"),
        phase_deviation=FloatSlider(
            value=np.pi, min=0.1, max=3 * np.pi, step=0.1,
            description="kp (rad)"),
    )


# ---------------------------------------------------------------------------
# Standalone execution
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    matplotlib.use("Agg")

    fig = plot_pm_demo()
    fig.savefig("pm_demo.png", dpi=150)
    print("Saved pm_demo.png")

    plt.show()
