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

# DSP & RF Fundamentals

A hands-on guide to Digital Signal Processing and Radio Frequency concepts.

**Audience**: Technical people (software devs, IT, etc.) who know nothing about RF or electronics.

Each section imports from standalone `.py` scripts that can also be run independently. Use the interactive sliders to explore how changing parameters affects the signals.

**Requirements**: `pip install numpy scipy matplotlib ipywidgets jupyter`


## Sections

Each section is a standalone notebook — open them in order or jump to whichever topic you want.

1. [Part 1 — Fundamentals](fundamentals/fundamentals.ipynb): signals, FFT, sampling, filtering
2. [Part 2 — Analog Modulation](modulations/analog_modulation.ipynb): AM, FM, PM
3. [Part 3 — Digital Modulation](modulations/digital_modulation.ipynb): ASK/PAM, BPSK, QPSK, FSK, QAM, BER
4. [Part 4 — Spread Spectrum](spreadspectrum/spread_spectrum.ipynb): LFSR, DSSS, FHSS, GPS Gold codes
5. [Part 5 — Antenna Arrays](arrays/antenna_arrays.ipynb): phased arrays, AOA



## Going deeper

The main sections stay intuition-first — no math prerequisites. If you want the theory behind a concept (complex numbers & I/Q, decibels, the FFT formula, BER curves), see the [**Deep Dive appendix**](deep_dive.ipynb). Individual sections also link to specific appendix entries via 🔍 *Going deeper* collapsibles.



---
# Summary

| Concept | Key Idea | Script |
|---------|----------|--------|
| Signals & FFT | Any signal = sum of sinusoids; FFT reveals them | `fundamentals/signal_basics.py` |
| Sampling | Must sample at 2x max frequency (Nyquist) | `fundamentals/sampling_aliasing.py` |
| Filtering | Keep desired frequencies, remove unwanted | `filters/filters.py` |
| Adding vs. multiplying | Sum stacks tones (DTMF); product creates sum & difference (mixing) | `fundamentals/mixing.py` |
| AM | Vary carrier amplitude to encode info | `modulations/am.py` |
| FM | Vary carrier frequency — more noise-resistant | `modulations/fm.py` |
| PM | Vary carrier phase — basis for digital PSK | `modulations/pm.py` |
| ASK/PAM | Digital amplitude modulation, eye diagrams | `modulations/ask_pam.py` |
| BPSK | 1 bit/symbol via phase flip — very robust | `modulations/bpsk.py` |
| QPSK | 2 bits/symbol via I/Q channels | `modulations/qpsk.py` |
| FSK | Frequency switching — robust to fading | `modulations/fsk.py` |
| QAM | Amplitude + phase — high throughput (WiFi, 5G) | `modulations/qam.py` |
| Channel & BER | Noise degrades signals; BER measures quality | `modulations/channel_ber.py` |
| LFSR | Pseudo-random sequences with autocorrelation | `spreadspectrum/lfsr_demo.py` |
| DSSS | Spread by multiplying with PN code | `spreadspectrum/dsss.py` |
| FHSS | Spread by hopping frequencies | `spreadspectrum/fhss.py` |
| GPS Gold Codes | CDMA — multiple satellites share one frequency | `spreadspectrum/gps-gold-codes/gps_gold_codes.py` |
| Phased Arrays | Multiple antennas for direction finding | `arrays/phased_array.py` |
