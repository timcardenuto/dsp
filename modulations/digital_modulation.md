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

# Part 3: Digital Modulation

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
# Part 3: Digital Modulation

Modern communications are digital — we transmit 1s and 0s, not continuous waveforms. Digital modulation maps discrete data (bits) onto a carrier wave.

The key tradeoff: **more bits per symbol = higher data rate but more susceptible to noise**.

| Modulation | Bits/Symbol | Robustness | Used In |
|------------|-------------|------------|----------|
| ASK/OOK | 1 | High | RFID, simple remotes |
| BPSK | 1 | Very High | GPS, deep space |
| QPSK | 2 | High | Satellite, early WiFi |
| FSK | 1 | High | Bluetooth, pagers |
| 16-QAM | 4 | Moderate | WiFi, 4G |
| 64-QAM | 6 | Low | WiFi, cable TV |
| 256-QAM | 8 | Very Low | WiFi 6, 5G |


## 3.1 ASK & PAM

**ASK** (Amplitude Shift Keying) is the digital version of AM — the simplest digital modulation. On-Off Keying (OOK) is the most basic form: carrier on = 1, carrier off = 0.

**PAM** (Pulse Amplitude Modulation) extends this to multiple amplitude levels. 4-PAM uses 4 levels to encode 2 bits per symbol.

The **eye diagram** is a powerful tool for visualizing signal quality — it overlays many symbol periods on top of each other. A wide-open "eye" means the signal is clean and easy to decode.

```python
from modulations.ask_pam import plot_ask_pam_demo
fig = plot_ask_pam_demo()
# fig.show()
```

```python
from modulations.ask_pam import interactive_ask_pam
interactive_ask_pam()
```

## 3.2 BPSK (Binary Phase Shift Keying)

BPSK flips the carrier's phase by 180 degrees to encode bits: 0 -> phase=0, 1 -> phase=180.

It's only 1 bit per symbol, but extremely robust against noise — used in GPS, deep-space communication, and any scenario where reliability matters more than speed.

```python
from modulations.bpsk import plot_bpsk_demo
fig = plot_bpsk_demo()
# fig.show()
```

```python
from modulations.bpsk import interactive_bpsk
interactive_bpsk()
```

## 3.3 QPSK (Quadrature Phase Shift Keying)

QPSK uses 4 phase states (45, 135, 225, 315 degrees) to encode **2 bits per symbol** — double the data rate of BPSK in the same bandwidth.

It works by splitting the carrier into two orthogonal channels:
- **I (In-phase)**: cosine component
- **Q (Quadrature)**: sine component

Each channel independently carries 1 bit, giving 2 bits per symbol total.


<details>
<summary>🔍 <b>Going deeper</b> — what "I/Q" and "quadrature" really mean (optional)</summary>

An I/Q sample is just an (x, y) coordinate of a rotating arrow (phasor). QPSK puts the arrow at 4 corners of a square; QAM uses a whole grid. See <a href="../deep_dive.ipynb#iq">Deep Dive → Complex numbers & phasors</a> if you want the full picture.

</details>


```python
from modulations.qpsk import plot_qpsk_demo
fig = plot_qpsk_demo()
# fig.show()
```

```python
from modulations.qpsk import interactive_qpsk
interactive_qpsk()
```

## 3.4 FSK (Frequency Shift Keying)

FSK is the digital version of FM — it switches between discrete frequencies to represent bits. Bit 0 uses one frequency, bit 1 uses another.

FSK is less spectrally efficient than PSK or QAM, but it's robust against amplitude fading (since information is in the frequency, not the amplitude). Used in Bluetooth, pagers, and low-power IoT.

```python
from modulations.fsk import plot_fsk_demo
fig = plot_fsk_demo()
# fig.show()
```

```python
from modulations.fsk import interactive_fsk
interactive_fsk()
```

## 3.5 QAM (Quadrature Amplitude Modulation)

QAM combines amplitude AND phase modulation — each symbol is a unique point on a 2D grid (the constellation diagram). 16-QAM uses a 4x4 grid to encode 4 bits per symbol; 64-QAM uses 8x8 for 6 bits.

QAM is the workhorse of modern communications: WiFi, 4G/5G, cable TV, and DSL all use it. Higher-order QAM packs more data into the same bandwidth, but requires a cleaner (higher SNR) channel.

```python
from modulations.qam import plot_qam_demo
fig = plot_qam_demo()
# fig.show()
```

```python
from modulations.qam import interactive_qam
interactive_qam()
```

## 3.6 Channel Effects & BER

In the real world, signals travel through a **channel** (air, cable, fiber) that adds noise and distortion. The key metric is **BER** (Bit Error Rate) — the fraction of bits received incorrectly.

- **SNR** (Signal-to-Noise Ratio): how much stronger the signal is vs the noise, in dB
- Higher SNR = fewer errors = lower BER
- Different modulations have different BER curves — BPSK is most robust, 64-QAM is most fragile

Watch the constellation degrade as SNR decreases — the clean grid of points becomes a noisy cloud.


<details>
<summary>🔍 <b>Going deeper</b> — SNR, Eb/N0, and reading BER curves (optional)</summary>

Higher-order modulations (16-QAM, 64-QAM) squeeze more bits per symbol but need more Eb/N0 to hit the same error rate. See <a href="../deep_dive.ipynb#ber">Deep Dive → SNR, Eb/N0, and the BER curve</a> for the math and a comparison plot.

</details>


```python
from modulations.channel_ber import plot_channel_demo
fig = plot_channel_demo()
# fig.show()
```

```python
from modulations.channel_ber import interactive_channel
interactive_channel()
```
