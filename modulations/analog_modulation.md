---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.19.4
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Part 2: Analog Modulation

Part of the [DSP & RF Fundamentals](../dsp_learning.ipynb) series.


```python
import sys, os
# Add parent dir (dsp/) so we can import fundamentals.*, modulations.*, etc.
_parent = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
if _parent not in sys.path:
    sys.path.insert(0, _parent)
%matplotlib inline

```

---
# Part 2: Analog Modulation

**Modulation** is the process of encoding information onto a high-frequency **carrier wave** so it can be transmitted through the air (or a cable).

Why not just transmit the original signal? Because low-frequency signals (like audio at 20-20,000 Hz) would need enormous antennas and would all interfere with each other. By shifting the signal up to a carrier frequency (e.g., 100 MHz for FM radio), we can use practical antennas and give each station its own slice of spectrum.

There are three properties of a sine wave we can vary to encode information:
- **Amplitude** → AM (Amplitude Modulation)
- **Frequency** → FM (Frequency Modulation)
- **Phase** → PM (Phase Modulation)


## 2.1 Amplitude Modulation (AM)

AM is the simplest modulation: the carrier's amplitude follows the message signal. It's how AM radio works — the audio signal controls the "volume" of the carrier wave.

The **modulation index** (0 to 1) controls how much the amplitude varies. If it exceeds 1.0 (overmodulation), the signal distorts.

```python
from modulations.am import plot_am_demo
fig = plot_am_demo()
# fig.show()
```

```python
from modulations.am import interactive_am
interactive_am()
```

### 🎧 Hear it

Plots show what modulation does; audio makes it *obvious*. Below: a clean 440 Hz tone, then what comes out after AM modulation + envelope demodulation at three modulation indices. At **m = 1.5** the envelope crosses zero and the demod distorts — you'll hear the buzz.


```python
from modulations._audio import audio, am_audio_demo
from IPython.display import display, Markdown

msg, results = am_audio_demo()
display(Markdown("**Original message (440 Hz tone):**"))
display(audio(msg))
for m, demod in results:
    label = "clean" if m < 1 else ("edge of over-modulation" if m == 1 else "over-modulated — hear the distortion")
    display(Markdown(f"**Demodulated at m = {m} ({label}):**"))
    display(audio(demod))

```

## 2.2 Frequency Modulation (FM)

FM encodes information by varying the carrier's *frequency* instead of its amplitude. When the message signal is positive, the carrier frequency increases; when negative, it decreases.

FM is more resistant to noise than AM because noise primarily affects amplitude, not frequency. This is why FM radio sounds cleaner than AM radio.

```python
from modulations.fm import plot_fm_demo
fig = plot_fm_demo()
# fig.show()
```

```python
from modulations.fm import interactive_fm
interactive_fm()
```

### 🎧 Hear it

FM's headline feature is **noise robustness**. Below is the same 440 Hz tone recovered from an FM signal at three SNR levels. AM would be unintelligible at 0 dB; FM is still audible.


```python
from modulations._audio import audio, fm_audio_demo
from IPython.display import display, Markdown

msg, results = fm_audio_demo()
display(Markdown("**Original message:**"))
display(audio(msg))
for snr, demod in results:
    display(Markdown(f"**FM demodulated at SNR = {snr} dB:**"))
    display(audio(demod))

```

## 2.3 Phase Modulation (PM)

PM encodes information by shifting the carrier's *phase* — the "starting position" of the wave cycle. PM and FM are closely related: PM of a signal is mathematically equivalent to FM of that signal's derivative.

PM is less common in analog systems but forms the basis for digital PSK (Phase Shift Keying), which is used in WiFi, cellular, and satellite communications.

```python
from modulations.pm import plot_pm_demo
fig = plot_pm_demo()
# fig.show()
```

```python
from modulations.pm import interactive_pm
interactive_pm()
```

### 🎧 Hear it

Increasing the phase deviation is the PM analog of turning up AM's modulation index — the recovered tone gets louder. Beyond π the demod starts to alias because phase wraps.


```python
from modulations._audio import audio, pm_audio_demo
from IPython.display import display, Markdown
import numpy as np

msg, results = pm_audio_demo()
display(Markdown("**Original message:**"))
display(audio(msg))
for pd, demod in results:
    display(Markdown(f"**PM demodulated, phase deviation = {pd/np.pi:.2f}π rad:**"))
    display(audio(demod))

```
