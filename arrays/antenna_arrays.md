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

# Part 5: Antenna Arrays

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
# Part 5: Antenna Arrays

A single antenna receives signals from all directions equally. A **phased array** uses multiple antennas to focus on a specific direction — like cupping your hand around your ear.

By measuring the phase differences between antennas, we can determine the **Angle of Arrival (AOA)** of an incoming signal. More antennas = narrower beam = better angular resolution.

```python
from arrays.phased_array import plot_array_demo
fig = plot_array_demo()
# fig.show()
```

```python
from arrays.phased_array import interactive_array
interactive_array()
```
