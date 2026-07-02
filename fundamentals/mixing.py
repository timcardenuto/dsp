"""Combining signals: adding vs. multiplying.

Adding two sinusoids just stacks them — you hear a chord and see two spikes in
the spectrum. Multiplying two sinusoids replaces the originals with their
**sum and difference frequencies** — the trigonometric identity

    cos(A) * cos(B) = 0.5 * [ cos(A - B) + cos(A + B) ]

is the seed of AM, RF mixers, and heterodyne receivers.

This module builds intuition with tones you can both see and hear.
"""

import numpy as np
import matplotlib.pyplot as plt

FS = 8000  # audio sample rate
DURATION = 1.5  # seconds — long enough to hear, short enough to plot


def tone(freq, duration=DURATION, fs=FS, amplitude=1.0):
    """Generate a pure sine tone."""
    t = np.arange(0, duration, 1 / fs)
    return t, amplitude * np.sin(2 * np.pi * freq * t)


def _plot_time_and_spectrum(t, sig, fs, title, xlim_ms=20, freq_max=None):
    """Two-panel figure: time (zoomed) and magnitude spectrum."""
    fig, (ax_t, ax_f) = plt.subplots(1, 2, figsize=(12, 3.2))
    fig.suptitle(title, fontsize=12)

    ax_t.plot(t * 1000, sig, linewidth=0.9)
    ax_t.set_xlim(0, xlim_ms)
    ax_t.set_xlabel("Time (ms)")
    ax_t.set_ylabel("Amplitude")
    ax_t.set_title("Time domain (first 20 ms)")
    ax_t.grid(True, alpha=0.3)

    N = len(sig)
    spectrum = np.abs(np.fft.rfft(sig)) / N * 2
    freqs = np.fft.rfftfreq(N, 1 / fs)
    ax_f.plot(freqs, spectrum, linewidth=0.9)
    ax_f.set_xlabel("Frequency (Hz)")
    ax_f.set_ylabel("Magnitude")
    ax_f.set_title("Frequency spectrum")
    ax_f.set_xlim(0, freq_max or fs / 2)
    ax_f.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


def demo_single(freq, fs=FS):
    """Return (fig, audio_array, fs) for a single tone."""
    t, sig = tone(freq, fs=fs)
    fig = _plot_time_and_spectrum(t, sig, fs, f"Single tone at {freq} Hz", freq_max=3000)
    return fig, sig, fs


def demo_sum(f1, f2, fs=FS):
    """Add two tones. Result: both spikes in the spectrum, chord in the audio."""
    t, s1 = tone(f1, fs=fs)
    _, s2 = tone(f2, fs=fs)
    summed = 0.5 * (s1 + s2)
    fig = _plot_time_and_spectrum(
        t, summed, fs,
        f"Adding {f1} Hz + {f2} Hz — spectrum has BOTH original frequencies",
        freq_max=3000,
    )
    return fig, summed, fs


def demo_product(f1, f2, fs=FS):
    """Multiply two tones. Result: the originals disappear; sum & difference appear.

    cos(2π f1 t) · cos(2π f2 t) = 0.5·[cos(2π(f1−f2)t) + cos(2π(f1+f2)t)]
    """
    t, s1 = tone(f1, fs=fs)
    _, s2 = tone(f2, fs=fs)
    product = s1 * s2
    fig = _plot_time_and_spectrum(
        t, product, fs,
        f"Multiplying {f1} Hz × {f2} Hz — spectrum shows |{f1}−{f2}| = {abs(f1-f2)} Hz "
        f"and {f1}+{f2} = {f1+f2} Hz (originals are GONE)",
        freq_max=3000,
    )
    return fig, product, fs


# ---------------------------------------------------------------------------
# DTMF — the canonical "two tones summed" system
# ---------------------------------------------------------------------------

DTMF_ROWS = {'1': 697, '2': 697, '3': 697,
             '4': 770, '5': 770, '6': 770,
             '7': 852, '8': 852, '9': 852,
             '*': 941, '0': 941, '#': 941}
DTMF_COLS = {'1': 1209, '2': 1336, '3': 1477,
             '4': 1209, '5': 1336, '6': 1477,
             '7': 1209, '8': 1336, '9': 1477,
             '*': 1209, '0': 1336, '#': 1477}


def dtmf_digit(digit, duration=0.25, fs=FS):
    """Generate one DTMF digit — the sum of its row and column frequencies."""
    f_row = DTMF_ROWS[digit]
    f_col = DTMF_COLS[digit]
    t = np.arange(0, duration, 1 / fs)
    return 0.5 * (np.sin(2 * np.pi * f_row * t) + np.sin(2 * np.pi * f_col * t))


def dtmf_sequence(digits, digit_duration=0.25, gap=0.08, fs=FS):
    """Concatenate DTMF tones with silence gaps — like dialing a phone."""
    gap_samples = np.zeros(int(gap * fs))
    parts = []
    for d in digits:
        parts.append(dtmf_digit(d, digit_duration, fs))
        parts.append(gap_samples)
    return np.concatenate(parts)


def demo_superhet(target_mhz=96.7, if_mhz=10.7,
                  stations=(88.5, 92.3, 96.7, 101.1, 105.9)):
    """Schematic of a superheterodyne receiver.

    Shows the antenna spectrum (several stations), the local oscillator, and
    the spectrum after multiplication — with an IF filter passband highlighting
    the station we tuned to.

    Not a numerical simulation — the frequencies are MHz and would need huge
    sample rates. This is a diagram of what happens *at each stage*, drawn
    from the sum/difference identity.
    """
    lo_mhz = target_mhz - if_mhz
    stations = np.asarray(stations, dtype=float)

    fig, axes = plt.subplots(3, 1, figsize=(12, 7))
    fig.suptitle(
        f"Superheterodyne receiver tuned to {target_mhz} MHz "
        f"(LO = {lo_mhz:.1f} MHz, IF = {if_mhz} MHz)",
        fontsize=12,
    )

    # ---------- Panel 1: antenna spectrum ----------
    ax = axes[0]
    for f in stations:
        color = "tab:red" if f == target_mhz else "tab:blue"
        lw = 2.5 if f == target_mhz else 1.2
        ax.stem([f], [1], linefmt=color, markerfmt=" ", basefmt=" ")
        ax.plot([f], [1], marker="v", color=color, markersize=8 if f == target_mhz else 5)
        ax.text(f, 1.08, f"{f}", ha="center", fontsize=8, color=color)
    ax.set_xlim(80, 115)
    ax.set_ylim(0, 1.35)
    ax.set_yticks([])
    ax.set_xlabel("Frequency (MHz)")
    ax.set_title("1. Antenna — every FM station arrives at once (target in red)")
    ax.grid(True, axis="x", alpha=0.3)

    # ---------- Panel 2: local oscillator ----------
    ax = axes[1]
    ax.stem([lo_mhz], [1], linefmt="tab:green", markerfmt=" ", basefmt=" ")
    ax.plot([lo_mhz], [1], marker="v", color="tab:green", markersize=8)
    ax.text(lo_mhz, 1.08, f"LO = {lo_mhz:.1f} MHz", ha="center", fontsize=9, color="tab:green")
    ax.set_xlim(80, 115)
    ax.set_ylim(0, 1.35)
    ax.set_yticks([])
    ax.set_xlabel("Frequency (MHz)")
    ax.set_title("2. Local oscillator — a single tone you control with the tuning dial")
    ax.grid(True, axis="x", alpha=0.3)

    # ---------- Panel 3: after multiplication ----------
    ax = axes[2]
    diff_freqs = np.abs(stations - lo_mhz)
    sum_freqs = stations + lo_mhz

    # Difference band (low side) — this is where the IF filter looks
    for f_station, f_diff in zip(stations, diff_freqs):
        color = "tab:red" if f_station == target_mhz else "tab:blue"
        lw = 2.5 if f_station == target_mhz else 1.2
        ax.stem([f_diff], [1], linefmt=color, markerfmt=" ", basefmt=" ")
        ax.text(f_diff, 1.05, f"{f_diff:.1f}", ha="center", fontsize=7, color=color)

    # Sum band (high side) — filtered out
    for f_station, f_sum in zip(stations, sum_freqs):
        ax.stem([f_sum], [0.6], linefmt="lightgray", markerfmt=" ", basefmt=" ")

    # IF filter passband
    passband_width = 0.4
    ax.axvspan(if_mhz - passband_width, if_mhz + passband_width,
               color="tab:green", alpha=0.2, label=f"IF filter @ {if_mhz} MHz")
    ax.axvline(if_mhz, color="tab:green", linestyle="--", alpha=0.6)

    ax.set_xlim(0, 220)
    ax.set_ylim(0, 1.35)
    ax.set_yticks([])
    ax.set_xlabel("Frequency (MHz)")
    ax.set_title(
        "3. After multiplying antenna × LO — each station shifts to |f − LO| and f + LO. "
        "Only the target lands in the IF filter."
    )
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, axis="x", alpha=0.3)

    plt.tight_layout()
    return fig


def demo_dtmf(digits="5", fs=FS):
    """Play a DTMF sequence with waveform + spectrogram.

    The spectrogram makes it visually clear that each key is two horizontal
    stripes — one row tone, one column tone.
    """
    sig = dtmf_sequence(digits, fs=fs)
    t = np.arange(len(sig)) / fs

    fig, (ax_t, ax_s) = plt.subplots(2, 1, figsize=(12, 5))
    fig.suptitle(f"DTMF for \"{digits}\" — each key = row tone + column tone (summed)",
                 fontsize=12)

    ax_t.plot(t, sig, linewidth=0.5)
    ax_t.set_xlabel("Time (s)")
    ax_t.set_ylabel("Amplitude")
    ax_t.set_title("Waveform")
    ax_t.grid(True, alpha=0.3)

    # Silence gaps between digits would produce exact zeros in the FFT
    # (and log10 warnings). Add inaudible broadband noise to avoid it.
    noise_floor = 1e-6 * np.random.default_rng(0).standard_normal(sig.shape)
    ax_s.specgram(sig + noise_floor, NFFT=512, Fs=fs, noverlap=384, cmap="magma")
    ax_s.set_ylim(500, 1700)
    ax_s.set_xlabel("Time (s)")
    ax_s.set_ylabel("Frequency (Hz)")
    ax_s.set_title("Spectrogram — two horizontal lines per digit")

    plt.tight_layout()
    return fig, sig, fs
