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

# Part 1: Fundamentals

Part of the [DSP & RF Fundamentals](../dsp_learning.ipynb) series.


```python
import sys, os
# Add parent dir (dsp/) so we can import fundamentals.*, modulations.*, etc.
_parent = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
if _parent not in sys.path:
    sys.path.insert(0, _parent)
%matplotlib inline
import matplotlib.pyplot as plt

```

---
# Part 1: Fundamentals

Before diving into radio, we need to understand three building blocks:
1. **Signals** — values that change over time (voltage, pressure, etc.)
2. **Frequency** — how fast a signal oscillates, measured in Hertz (cycles/second)
3. **The FFT** — an algorithm that decomposes any signal into its constituent frequencies

Think of the FFT like a prism splitting white light into a rainbow — it reveals the hidden frequency components inside a complex signal.


## 1.1 Signal Basics & FFT

Any signal can be built by adding sinusoids (sine waves) together. Each sinusoid has:
- **Frequency**: how fast it oscillates
- **Amplitude**: how tall the peaks are
- **Phase**: where in its cycle it starts

The plots below show three sinusoids combined into one composite signal, then the FFT reveals which frequencies are present.


<details>
<summary>🔍 <b>Going deeper</b> — what the FFT is actually computing (optional)</summary>

The FFT is a correlation: for each frequency bin, it checks how well a reference sinusoid at that frequency matches your signal. If it matches, the sum is big. See <a href="../deep_dive.ipynb#fft">Deep Dive → The FFT in one page</a> for the formula and what "windowing" means.

</details>


```python
from fundamentals.signal_basics import plot_signal_basics, plot_spectrogram_demo
fig = plot_signal_basics()
# fig.show()
```

**Try it yourself** — drag the sliders to change frequencies and amplitudes. Watch how the time-domain waveform and frequency spectrum change together.

```python
from fundamentals.signal_basics import interactive_signal
interactive_signal()
```

### Spectrogram

A spectrogram shows how frequencies change *over time* — it's a 2D heatmap with time on the x-axis, frequency on the y-axis, and color representing power. Below is a "chirp" signal that sweeps from 200 Hz to 2000 Hz.

```python
fig = plot_spectrogram_demo()
# fig.show()
```


## 1.2 Combining signals: adding vs. multiplying

Two of the most important operations in DSP look almost identical on the page but do wildly different things to a signal.

- **Adding two signals** stacks them. You *hear both at once* (a chord), and the spectrum shows a spike at each original frequency. This is how DTMF (touch-tone phones), musical chords, and OFDM subcarriers work.
- **Multiplying two signals** is called **mixing**. The originals *disappear* and get replaced by their **sum** and **difference** frequencies. This is the mathematical seed of AM modulation, RF mixers, and every heterodyne receiver.

The reason multiplication has that effect is one trig identity:

$$\cos(A) \cdot \cos(B) = \tfrac{1}{2}\bigl[\cos(A - B) + \cos(A + B)\bigr]$$

Below: see it and hear it.



### First, hear each tone by itself

A pure 440 Hz sine (musical A) and a 660 Hz sine (roughly E, a fifth above).

```python
from fundamentals.mixing import demo_single
from IPython.display import display, Audio, Markdown

for freq in (440, 660):
    fig, sig, fs = demo_single(freq)
    display(Markdown(f"**{freq} Hz — single tone**"))
    display(fig); plt.close(fig)
    display(Audio(sig, rate=fs, normalize=True))

```

### Adding them — you hear both at once

Notice: the time-domain waveform looks messy but the spectrum shows **two clean spikes at the original frequencies**. Nothing new was created. Your ear hears a two-note chord.

```python
from fundamentals.mixing import demo_sum

fig, sig, fs = demo_sum(440, 660)
display(fig); plt.close(fig)
display(Audio(sig, rate=fs, normalize=True))

```

### Multiplying them — the frequencies shift

Watch the spectrum carefully. The 440 Hz and 660 Hz spikes are **gone**. In their place: a spike at 660 − 440 = **220 Hz** and one at 660 + 440 = **1100 Hz**. That's the mixing identity in action.

**One of the most important tricks in RF.** Multiplication *moves signal energy to new frequencies* — specifically to the sum and difference. Every AM transmitter, every RF mixer, every superheterodyne receiver relies on this. Want to shift a 1 kHz audio signal up to 100 MHz for transmission? Multiply it by a 100 MHz carrier — energy appears at 100 MHz ± 1 kHz. Want to shift a received 2.4 GHz signal down to something your ADC can sample? Multiply by a 2.399 GHz local oscillator — the signal appears at 1 MHz (the difference) and 4.799 GHz (the sum, which you filter out).

The audio isn't the point here — the frequency shift is. But if you're curious why our two examples *sound* so different, see the footnote below.

> Forward pointer: AM modulation is literally `message × carrier`. See [Part 2: Analog Modulation](../modulations/analog_modulation.ipynb).

---

<details>
<summary>📝 <b>Footnote</b> — why the sum sounds like a chord and the product sounds like something else</summary>

Both cases produce **two sinusoids playing together**. Your ear does its own spectral analysis and hears whatever frequencies are present — there's no special "multiplication timbre." The perceived difference is entirely due to *which* frequencies ended up in each output:

- **Sum (440 + 660 Hz)**: a 3:2 frequency ratio — a **perfect fifth**, one of the most consonant intervals in music. Both notes sit in the musical mid-range. Your ear/brain fuses them into a chord.
- **Product (220 + 1100 Hz)**: a 5:1 ratio, widely spaced and dissonant. The 220 Hz gives a bass buzz; the 1100 Hz sits in the harsh sibilance range. Wide, non-integer spacing is exactly what makes bells and gongs sound bell-like — their partials aren't harmonically related either.

So the "ring modulator" effect people associate with multiplication is really about the *choice of inputs*: multiply a harmonically rich signal (voice, guitar) by a sine and every harmonic gets shifted, creating a whole new inharmonic spectrum. With two pure sines, you just get two pure sines out — a boring truth the marketing of "ring mod" tends to obscure.

</details>


```python
from fundamentals.mixing import demo_product

fig, sig, fs = demo_product(440, 660)
display(fig); plt.close(fig)
display(Audio(sig, rate=fs, normalize=True))

```

### Dual-Tone Multi-Frequency (DTMF) — addition in the real world

When you press "5" on a phone, it plays sin(770 Hz) + sin(1336 Hz). No multiplication, just addition — two tones summed. The spectrogram below shows exactly two horizontal stripes per digit, one for the row frequency and one for the column frequency.

Try changing `digits` to a phone number and re-running.

```python
from fundamentals.mixing import demo_dtmf

digits = "5551234"
fig, sig, fs = demo_dtmf(digits)
display(fig); plt.close(fig)
display(Audio(sig, rate=fs, normalize=True))

```

## 1.3 Filtering

Filters let you keep the frequencies you want and remove the ones you don't:
- **Low-pass**: keeps low frequencies, blocks high (like bass boost)
- **High-pass**: keeps high frequencies, blocks low (like treble boost)
- **Band-pass**: keeps a range of frequencies, blocks everything else (like tuning a radio to one station)

Filtering is an important operation in DSP — it's used everywhere from noise removal to signal recovery.

```python
from filters.filters import plot_filter_demo
fig = plot_filter_demo()
# fig.show()
```

```python
from filters.filters import interactive_filter
interactive_filter()
```

## 1.4 Sampling & Aliasing

Real-world signals are continuous, but computers work with discrete samples. **Sampling** is the process of measuring a signal at regular intervals.

The **Nyquist theorem** says: to capture a signal faithfully, you must sample at least **twice** its highest frequency. If you sample too slowly, the signal "folds" to the wrong frequency — this is **aliasing**.

Think of aliasing like a car wheel in a movie — if the camera frame rate is too slow relative to the wheel's rotation, the wheel appears to spin backward.

```python
from fundamentals.sampling_aliasing import plot_sampling_demo
fig = plot_sampling_demo()
# fig.show()
```

```python
from fundamentals.sampling_aliasing import interactive_sampling
interactive_sampling()
```



### Superheterodyne — multiplication in the wild

DTMF was addition in the real world. Multiplication's equivalent is the **superheterodyne receiver**, the design in every FM radio, phone, TV, and GPS receiver since Edwin Armstrong invented it in 1918.

**The problem:** your antenna picks up every FM station simultaneously (88–108 MHz), plus a lot of noise. You want to hear *one* station. Building a filter sharp enough to pick out just 96.7 MHz while rejecting 96.5 and 96.9 — and also being able to *retune* it — is extremely hard.

**The trick:** don't move the filter, move the signal. Build one very good filter at a fixed **intermediate frequency (IF)** — 10.7 MHz for FM. Then use multiplication to shift whichever station you want down to that filter.

**How:** the receiver contains a **local oscillator (LO)** whose frequency you control with the tuning dial. It multiplies the antenna signal by the LO. From the sum/difference rule, each station at frequency $f$ shows up at $|f - \text{LO}|$ and $f + \text{LO}$. Set the LO so `target − LO = 10.7 MHz` and the target lands *exactly* in the IF filter. Every other station lands elsewhere and gets filtered out.

**Turning the dial changes the LO frequency.** That's it. That's what tuning a radio does.

Below: a schematic of the three stages when tuned to 96.7 MHz. Try changing `target` and re-running to "turn the dial".

```python
from fundamentals.mixing import demo_superhet

target = 96.7   # MHz — change this to tune to a different station
fig = demo_superhet(target_mhz=target)
display(fig); plt.close(fig)

```
