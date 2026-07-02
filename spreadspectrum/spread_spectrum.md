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

# Part 4: Spread Spectrum

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
# Part 4: Spread Spectrum

Spread spectrum techniques deliberately spread a signal across a much wider bandwidth than necessary. This seems wasteful, but it provides:
- **Jam resistance**: a jammer must cover the entire wide band, not just one frequency
- **Signal hiding**: spread signals can be below the noise floor and still be recovered
- **Multiple access**: multiple users can share the same band simultaneously (CDMA)

Two main approaches:
- **DSSS** (Direct Sequence): multiply by a fast pseudo-random code to spread continuously
- **FHSS** (Frequency Hopping): rapidly switch between narrow frequencies


## 4.1 LFSR & Pseudo-Random Sequences

Before we can spread a signal, we need pseudo-random sequences. An **LFSR** (Linear Feedback Shift Register) generates these by shifting bits through a register and XOR-ing selected "tap" positions to produce the next input bit.

With the right taps, an n-bit LFSR produces a **maximal-length sequence** of 2^n - 1 bits before repeating. These sequences have a critical property: their **autocorrelation** shows a sharp peak at zero shift and is nearly flat everywhere else — essential for distinguishing signals in spread spectrum.

```python
from spreadspectrum.lfsr_demo import plot_lfsr_demo
fig = plot_lfsr_demo()
# fig.show()
```

```python
from spreadspectrum.lfsr_demo import interactive_lfsr
interactive_lfsr()
```

## 4.2 DSSS (Direct Sequence Spread Spectrum)

DSSS multiplies your data by a fast pseudo-random (PN) sequence. Each data bit is multiplied by many PN "chips", spreading it across a wider bandwidth.

The receiver multiplies by the same PN sequence to **despread** — the original signal pops back out, while interference from other sources stays spread (looks like low-level noise).

**Processing gain** = chips per bit. More chips = wider spread = more robust, but uses more bandwidth.

```python
from spreadspectrum.dsss import plot_dsss_demo
fig = plot_dsss_demo()
# fig.show()
```

```python
from spreadspectrum.dsss import interactive_dsss
interactive_dsss()
```

## 4.3 FHSS (Frequency Hopping Spread Spectrum)

Instead of spreading across all frequencies simultaneously (DSSS), FHSS rapidly **hops** between different narrow frequencies in a pseudo-random pattern.

The transmitter and receiver share the same hopping pattern. A jammer or eavesdropper who doesn't know the pattern can only interfere with a fraction of the hops.

Bluetooth uses FHSS — it hops between 79 channels, 1600 times per second.

```python
from spreadspectrum.fhss import plot_fhss_demo
fig = plot_fhss_demo()
# fig.show()
```

```python
from spreadspectrum.fhss import interactive_fhss
interactive_fhss()
```

## 4.4 GPS & Gold Codes

GPS is DSSS in action. Each satellite broadcasts on the **same frequency** but uses a unique 1023-chip **Gold code** (a special type of PN sequence with excellent cross-correlation properties).

Your GPS receiver knows all the codes. To find satellite #5, it correlates the received signal against PRN #5's code — that satellite's signal pops out, while all other satellites remain as low-level noise. This is **CDMA** (Code Division Multiple Access).

```python
sys.path.insert(0, os.path.join(os.getcwd(), 'gps-gold-codes'))
from gps_gold_codes import plot_gps_demo, plot_correlation_demo

fig1, fig2 = plot_gps_demo()
# fig1.show()
# fig2.show()

```

```python
fig = plot_correlation_demo()
# fig.show()
```

```python
from gps_gold_codes import interactive_gps
interactive_gps()
```
