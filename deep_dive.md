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

# Deep Dive Appendix

The main notebooks are intuition-first — no math prerequisites. This appendix is where the math lives, for when you want to understand *why* the intuition works.

You don't need to read this front-to-back. Sections are self-contained; the main notebooks link here from `Going deeper` collapsibles when a concept is relevant.

**Contents**
- [Complex numbers & phasors (what I/Q really is)](#iq)
- [Decibels (dB) — the log-ratio cheat sheet](#db)
- [The FFT in one page](#fft)
- [SNR, Eb/N0, and the BER curve](#ber)


```python
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```

<a id="iq"></a>
## Complex numbers & phasors — what I/Q really is

Every sine wave has three properties: **amplitude**, **frequency**, and **phase**. Tracking all three with plain trigonometry gets messy fast. Engineers instead picture a sine wave as a **rotating arrow** in a 2D plane:

- Length of the arrow = amplitude.
- How fast it spins = frequency.
- Starting angle = phase.

That arrow is called a **phasor**. Its horizontal coordinate is called **I** (in-phase); its vertical coordinate is **Q** (quadrature). When you project the tip onto the horizontal axis and watch that projection over time, you get a cosine wave. Onto the vertical axis, a sine wave.

Complex numbers are just a compact way to write down that (I, Q) pair: `I + jQ`. The identity that ties it all together is Euler's formula:

$$e^{j\omega t} = \cos(\omega t) + j\sin(\omega t)$$

Left side: "an arrow of length 1 spinning at rate ω." Right side: "its I coordinate is cos, its Q coordinate is sin." Same object, two notations.

**Why this matters practically:**
- **QPSK, QAM, and every I/Q-based modulation** work by placing a phasor at specific (I, Q) points to encode bits. QPSK uses 4 corners of a square; 16-QAM uses a 4×4 grid.
- **Software Defined Radios** output samples as I/Q pairs. Each sample *is* one snapshot of the phasor.
- **Multiplication of complex exponentials adds their angles**, which is why frequency mixing (heterodyning) is so clean: multiply by $e^{-j\omega_c t}$ and the carrier disappears, leaving just the baseband.

```python
# Visualize: rotating phasor and the sine wave it traces out.
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))

theta = np.linspace(0, 2 * np.pi, 200)
ax1.plot(np.cos(theta), np.sin(theta), '--', color='lightgray')
for angle in np.linspace(0, 2 * np.pi, 8, endpoint=False):
    ax1.arrow(0, 0, 0.9 * np.cos(angle), 0.9 * np.sin(angle),
              head_width=0.05, alpha=0.35, color='tab:blue', length_includes_head=True)
ax1.arrow(0, 0, 0.9 * np.cos(np.pi/6), 0.9 * np.sin(np.pi/6),
          head_width=0.08, color='tab:red', length_includes_head=True, linewidth=2)
ax1.set_xlim(-1.3, 1.3); ax1.set_ylim(-1.3, 1.3); ax1.set_aspect('equal')
ax1.axhline(0, color='k', linewidth=0.5); ax1.axvline(0, color='k', linewidth=0.5)
ax1.set_xlabel('I (real / cos)'); ax1.set_ylabel('Q (imag / sin)')
ax1.set_title('Phasor rotating at ω (red = one snapshot)')

t = np.linspace(0, 2, 500)
ax2.plot(t, np.cos(2 * np.pi * t), label='I = cos(ωt)')
ax2.plot(t, np.sin(2 * np.pi * t), label='Q = sin(ωt)')
ax2.set_xlabel('time'); ax2.set_title('Projections onto I and Q axes over time')
ax2.legend(); ax2.grid(True, alpha=0.3)
plt.tight_layout(); plt.show()
```

<a id="db"></a>
## Decibels (dB) — the log-ratio cheat sheet

RF and DSP use decibels everywhere. A dB is just a compressed way to write a ratio using logarithms — it turns huge multiplicative ranges (a signal a million times weaker than another) into small additive numbers.

For a **power** ratio:

$$X_{\text{dB}} = 10 \cdot \log_{10}\left(\frac{P_1}{P_2}\right)$$

For an **amplitude / voltage** ratio (power ∝ amplitude², so the coefficient becomes 20):

$$X_{\text{dB}} = 20 \cdot \log_{10}\left(\frac{V_1}{V_2}\right)$$

**Landmarks worth memorizing:**

| dB | Power ratio | Amplitude ratio | Feel |
|----|-------------|-----------------|------|
| 0 | 1× | 1× | Same |
| 3 | 2× | 1.41× | Doubled power |
| 6 | 4× | 2× | Doubled amplitude |
| 10 | 10× | 3.16× | One order of magnitude |
| 20 | 100× | 10× | Two orders |
| 30 | 1000× | 31.6× | Three orders |
| -3 | 0.5× | 0.71× | Half power (the classic filter cutoff) |

**Why dB is used, not raw ratios:**
- Filter cutoffs, antenna gains, and noise floors span 6+ orders of magnitude. On a linear scale you can't plot a receiver's dynamic range on one graph; on dB, it fits.
- Gains in cascaded stages *add* in dB instead of multiplying. A 10 dB amp followed by a 6 dB amp gives 16 dB total gain — no calculator required.

**dBm** is a specific case: power referenced to 1 milliwatt. `0 dBm = 1 mW`, `-30 dBm = 1 µW`, `+30 dBm = 1 W`.


<a id="fft"></a>
## The FFT in one page

The intuition given in Part 1 — "the FFT decomposes any signal into sinusoids" — is enough to use it. Here's what's actually happening:

For a sampled signal $x[n]$ with $N$ samples, the Discrete Fourier Transform (DFT) is:

$$X[k] = \sum_{n=0}^{N-1} x[n] \cdot e^{-j 2\pi k n / N}$$

Read this as: **for each frequency bin `k`, multiply the input by a reference sinusoid at that frequency and add up the products.** If the input contains that frequency, the products align in phase and the sum is large. If it doesn't, they cancel and the sum is near zero. It's a correlation, done once per bin.

The **F**ast in FFT is an algorithmic trick that computes the same result in $O(N \log N)$ instead of $O(N^2)$ by reusing subexpressions. Same output, faster.

**Practical things to know:**
- Bin `k` corresponds to frequency $k \cdot f_s / N$ Hz. Bin resolution improves with either lower sample rate or more samples.
- The output is complex. `|X[k]|` is the magnitude (how much of that frequency is present); `angle(X[k])` is the phase.
- For real input, the output is symmetric — bin `N-k` is the complex conjugate of bin `k`. Only bins `0` through `N/2` carry unique info. That's why real-signal spectrum plots show only the first half.
- **Windowing** matters. A raw FFT assumes the signal is periodic over your capture. Multiplying by a window function (Hann, Hamming, Blackman) before the FFT smooths the ends and reduces spectral leakage.



<a id="ber"></a>
## SNR, Eb/N0, and the BER curve

**SNR** is signal-to-noise ratio, usually expressed in dB. Higher is better; when SNR is low, bit errors happen.

**Eb/N0** ("ebb-no") is a per-bit version: energy per bit divided by noise power spectral density. It's the fair way to compare modulations — SNR is deceptive because a wideband modulation looks noisier than a narrowband one at the same information rate.

$$E_b/N_0 = \text{SNR} + 10 \log_{10}\!\left(\frac{B}{R_b}\right)$$

where $B$ is the receiver bandwidth and $R_b$ is the bit rate.

The **BER curve** plots bit error rate vs. Eb/N0 for a given modulation. Key patterns:
- Curves fall off exponentially — a couple of extra dB can swing BER by orders of magnitude. That's the "waterfall" shape.
- **More bits per symbol → worse BER at the same Eb/N0.** QPSK ties BPSK on Eb/N0 (a lucky geometry), but 16-QAM needs ~4 dB more, 64-QAM another ~5 dB, etc. This is the fundamental throughput-vs-robustness tradeoff.
- **Shannon's limit** ($E_b/N_0 = -1.59$ dB) is the theoretical floor. Practical systems approach it with coding (LDPC, turbo codes) but never cross it.

For the coherent BPSK / QPSK curve specifically:
$$P_b = Q\!\left(\sqrt{2 E_b / N_0}\right)$$
where $Q$ is the Q-function (tail of the standard normal). You don't need to derive this — just know that it's the reference every other modulation is compared against.

```python
# BER waterfall for common modulations
from scipy.special import erfc

def q(x):
    return 0.5 * erfc(x / np.sqrt(2))

ebno_db = np.linspace(0, 20, 200)
ebno = 10 ** (ebno_db / 10)

curves = {
    'BPSK / QPSK': q(np.sqrt(2 * ebno)),
    '16-QAM':      (3 / 4) * q(np.sqrt((4 / 5) * ebno)),
    '64-QAM':      (7 / 12) * q(np.sqrt((2 / 7) * ebno)),
}

fig, ax = plt.subplots(figsize=(9, 5))
for label, ber in curves.items():
    ax.semilogy(ebno_db, ber, label=label)
ax.set_xlabel('Eb/N0 (dB)'); ax.set_ylabel('Bit Error Rate')
ax.set_title('BER waterfall — more bits/symbol needs more Eb/N0')
ax.set_ylim(1e-6, 1)
ax.grid(True, which='both', alpha=0.3); ax.legend()
plt.show()
```
