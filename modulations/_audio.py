"""Audio helpers for the analog-modulation notebook.

Lets the reader *hear* what modulation and demodulation do — over-modulation
distortion, FM noise robustness, etc. — instead of only seeing them.

Everything uses audio-band sample rates (fs = 8000 Hz) and audio-band
frequencies so the resulting waveforms are playable through a browser via
IPython.display.Audio.
"""

import numpy as np
from IPython.display import Audio

from .am import am_modulate, am_demodulate
from .fm import fm_modulate, fm_demodulate
from .pm import pm_modulate, pm_demodulate


FS = 8000  # audio sample rate — small enough to keep vectors short


def audio(sig, fs=FS, normalize=True):
    """Return an Audio widget for `sig`. Normalizes to avoid clipping."""
    sig = np.asarray(sig, dtype=float)
    if normalize:
        peak = np.max(np.abs(sig))
        if peak > 0:
            sig = sig / peak * 0.9
    return Audio(sig, rate=fs, normalize=False)


def _message(fs=FS, duration=1.5, freq=440.0):
    """A clean 440 Hz tone (A above middle C) — the 'audio' we want to send."""
    t = np.arange(0, duration, 1 / fs)
    return t, np.sin(2 * np.pi * freq * t)


def am_audio_demo(mod_indices=(0.5, 1.0, 1.5), carrier_freq=2000, fs=FS):
    """Return (message, [(m, demodulated) for each mod index]).

    At m=0.5 the demod sounds clean, m=1.0 is the edge, m=1.5 is
    over-modulated — you can *hear* the buzz from envelope zero-crossings.
    """
    t, msg = _message(fs=fs)
    results = []
    for m in mod_indices:
        _, modulated, _ = am_modulate(msg, carrier_freq, fs, mod_index=m)
        demod = am_demodulate(modulated, fs, carrier_freq)
        # Remove DC so we hear the tone, not a click.
        demod = demod - np.mean(demod)
        results.append((m, demod))
    return msg, results


def fm_audio_demo(snr_dbs=(30, 10, 0), carrier_freq=2000, freq_deviation=500, fs=FS):
    """Return (clean_message, [(snr_db, demodulated) ...]).

    FM's famous noise robustness: at 30 dB the tone is crystal clear; at
    0 dB AM would be gone but FM is still audible.
    """
    t, msg = _message(fs=fs)
    fm_sig = fm_modulate(msg, carrier_freq, fs, freq_deviation)
    results = []
    rng = np.random.default_rng(0)
    for snr_db in snr_dbs:
        sig_power = np.mean(fm_sig ** 2)
        noise_power = sig_power / (10 ** (snr_db / 10))
        noisy = fm_sig + rng.normal(0, np.sqrt(noise_power), size=fm_sig.shape)
        demod = fm_demodulate(noisy, fs)
        demod = demod - np.mean(demod)
        results.append((snr_db, demod))
    return msg, results


def pm_audio_demo(phase_devs=(np.pi / 4, np.pi / 2, np.pi), carrier_freq=2000, fs=FS):
    """Return (message, [(phase_dev, demodulated) ...]).

    Small phase deviation sounds quiet, larger deviation sounds louder —
    directly analogous to raising the modulation index in AM but by
    wiggling phase instead of amplitude.
    """
    t, msg = _message(fs=fs)
    results = []
    for pd in phase_devs:
        _, pm_sig = pm_modulate(msg, carrier_freq, fs, pd)
        _, demod = pm_demodulate(pm_sig, fs, carrier_freq)
        demod = demod - np.mean(demod)
        results.append((pd, demod))
    return msg, results
