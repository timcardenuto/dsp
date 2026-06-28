"""
Linear Feedback Shift Registers (LFSRs) & Correlation Properties

An LFSR is a shift register whose input bit is computed by XOR-ing selected
bits ("taps") of the current register state.  Despite being completely
deterministic, the output sequence *looks* random -- it passes many
statistical tests for randomness.  This makes LFSRs cheap, fast sources of
pseudo-random bits used everywhere from GPS to Wi-Fi to scrambling pay-TV.

How it works (no hardware background needed):
  1. Start with a row of N bits (the "register") loaded with some non-zero
     initial value (the "seed").
  2. Each clock tick:
       a. Pick certain bit positions (the "taps") and XOR them together to
          get one new bit.
       b. Shift every bit one position to the right; the new bit enters on
          the left.
       c. The bit that falls off the right end is the output.
  3. Repeat.  The output is one bit per clock tick.

The magic comes from choosing the *right* taps.  If the tap positions
correspond to a "primitive polynomial," the register cycles through every
possible non-zero state before repeating.  For an N-bit register that means
2^N - 1 states -- the longest possible cycle.  These are called
"maximal-length sequences" (or "m-sequences").

Why we care -- correlation:
  - Autocorrelation: shift an m-sequence against itself.  At zero shift you
    get a sharp spike (every bit matches); at any other shift the correlation
    drops to -1/N, practically zero.  This "thumbtack" autocorrelation is
    what lets a GPS receiver lock onto a satellite signal buried in noise.
  - Cross-correlation: compare two *different* m-sequences.  The values stay
    low for all shifts, meaning the sequences are nearly orthogonal.  This
    lets multiple transmitters share the same frequency band (CDMA).

Polynomial notation (by example):
  "1 + x + x^4"  means tap positions 0 and 1 feed back into a 4-bit
  register.  The highest power (4) gives the register length; the lower
  powers tell you which bits to XOR together for the feedback bit.

This module provides clean functions for generating LFSR sequences and
computing their correlation, plus visualizations that show why these
properties matter for spread-spectrum communications.
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Core functions (importable by notebooks)
# ---------------------------------------------------------------------------

def lfsr_sequence(taps, seed, length):
    """Generate a pseudo-random bit sequence from a Linear Feedback Shift Register.

    The register is clocked *length* times.  Each tick the feedback bit is
    computed by XOR-ing the tap positions, then the register shifts right and
    the feedback bit enters on the left.  The output bit is the rightmost bit
    that shifts out.

    Parameters
    ----------
    taps : list of int
        Register positions (0-indexed) that are XOR-ed with the output bit
        to produce the feedback.  The output bit (rightmost register cell)
        is always included implicitly -- do NOT list it.  For example, the
        polynomial 1 + x + x^4 has one intermediate term x^1; the
        corresponding register position is 0, so pass ``taps=[0]``.
    seed : list or ndarray of int
        Initial register state, e.g. ``[1, 0, 0, 1]``.  Must be non-zero.
        The length of *seed* determines the register size.
    length : int
        Number of output bits to generate.  For one full maximal-length
        period, use ``2**len(seed) - 1``.

    Returns
    -------
    sequence : ndarray of int (0s and 1s), shape (length,)
    """
    order = len(seed)
    reg = np.array(seed, dtype=int)
    seq = np.empty(length, dtype=int)

    for i in range(length):
        # Output: rightmost bit
        seq[i] = reg[-1]

        # Feedback: output bit XOR'd with each tap position
        feedback = reg[-1]
        for t in taps:
            feedback ^= reg[t]

        # Shift right; feedback enters on the left
        reg = np.roll(reg, 1)
        reg[0] = feedback

    return seq


def autocorrelation(sequence):
    """Compute the circular (periodic) autocorrelation of a bipolar sequence.

    The input should already be in bipolar form (-1/+1).  The output is
    normalized so the peak at zero lag equals 1.0.

    Parameters
    ----------
    sequence : array-like of -1/+1 values

    Returns
    -------
    lags : ndarray – shift amounts (0 to N-1)
    corr : ndarray – normalized correlation at each lag
    """
    seq = np.asarray(sequence, dtype=float)
    N = len(seq)
    corr = np.empty(N)
    for k in range(N):
        corr[k] = np.dot(seq, np.roll(seq, k))
    corr = corr / N  # normalize so peak = 1.0
    lags = np.arange(N)
    return lags, corr


def cross_correlation(seq1, seq2):
    """Compute the circular cross-correlation of two bipolar sequences.

    Both inputs should be bipolar (-1/+1) and the same length.

    Parameters
    ----------
    seq1, seq2 : array-like of -1/+1 values

    Returns
    -------
    lags : ndarray – shift amounts (0 to N-1)
    corr : ndarray – normalized cross-correlation at each lag
    """
    s1 = np.asarray(seq1, dtype=float)
    s2 = np.asarray(seq2, dtype=float)
    N = len(s1)
    corr = np.empty(N)
    for k in range(N):
        corr[k] = np.dot(s1, np.roll(s2, k))
    corr = corr / N
    lags = np.arange(N)
    return lags, corr


def to_bipolar(sequence):
    """Map a unipolar (0/1) sequence to bipolar (-1/+1).

    0 -> -1,  1 -> +1.  This is the standard mapping used before computing
    correlation so that matched bits contribute +1 and mismatched bits -1.
    """
    return 2 * np.asarray(sequence, dtype=float) - 1


# ---------------------------------------------------------------------------
# Well-known maximal-length tap sets (polynomial -> tap positions)
# ---------------------------------------------------------------------------

# Each entry: order -> list of tap indices (0-indexed, excluding the order
# itself).  These correspond to primitive polynomials and produce sequences
# of length 2^order - 1.
MAXIMAL_LENGTH_TAPS = {
    3: [0],          # 1 + x + x^3
    4: [0],          # 1 + x + x^4
    5: [1],          # 1 + x^2 + x^5
    6: [0],          # 1 + x + x^6
    7: [2],          # 1 + x^3 + x^7
}

# Alternate taps that also produce maximal-length sequences (useful for
# generating a *different* sequence of the same length for cross-correlation
# demos).
ALTERNATE_TAPS = {
    4: [2],          # 1 + x^3 + x^4  (the "reversed" polynomial)
    5: [2],          # 1 + x^3 + x^5
}


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_lfsr_demo():
    """Generate a multi-panel figure demonstrating LFSR properties.

    Panel 1: Raw LFSR output (bipolar stem plot) for a 4th-order sequence.
    Panel 2: Autocorrelation of a maximal-length sequence -- sharp peak at
             zero shift, near-zero everywhere else.
    Panel 3: Cross-correlation between two different 4th-order m-sequences
             -- low and flat, showing near-orthogonality.
    Panel 4: Comparison of 3rd, 4th, and 5th order sequence lengths,
             illustrating exponential growth of the period.
    """
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle("Linear Feedback Shift Registers (LFSRs)", fontsize=14)

    # -- Panel 1: LFSR output sequence (4th order, 1 + x + x^4) ----------
    order = 4
    N = 2 ** order - 1
    seed = [1, 0, 0, 0]
    seq_bits = lfsr_sequence(MAXIMAL_LENGTH_TAPS[order], seed, N)
    seq_bp = to_bipolar(seq_bits)

    ax = axes[0, 0]
    ax.stem(np.arange(N), seq_bp, linefmt="tab:blue", markerfmt="o",
            basefmt="k-")
    ax.set_title("4th-Order LFSR Output (bipolar)")
    ax.set_xlabel("Chip index")
    ax.set_ylabel("Value (-1 / +1)")
    ax.set_ylim(-1.5, 1.5)
    ax.set_yticks([-1, 0, 1])
    ax.grid(True, alpha=0.3)
    ax.annotate(f"Period = 2^{order} - 1 = {N} chips",
                xy=(0.02, 0.92), xycoords="axes fraction", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow",
                          alpha=0.8))

    # -- Panel 2: Autocorrelation ----------------------------------------
    lags, acorr = autocorrelation(seq_bp)

    ax = axes[0, 1]
    ax.stem(lags, acorr, linefmt="tab:green", markerfmt="o", basefmt="k-")
    ax.set_title("Autocorrelation of 4th-Order m-Sequence")
    ax.set_xlabel("Circular shift (chips)")
    ax.set_ylabel("Normalized correlation")
    ax.grid(True, alpha=0.3)
    ax.axhline(y=-1/N, color="tab:red", linestyle="--", linewidth=0.8,
               label=f"-1/N = {-1/N:.3f}")
    ax.legend(fontsize=8)
    ax.annotate("Sharp peak at zero shift\n(key property for GPS/CDMA)",
                xy=(0.35, 0.75), xycoords="axes fraction", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow",
                          alpha=0.8))

    # -- Panel 3: Cross-correlation --------------------------------------
    # Generate a second 4th-order sequence with different taps
    seq2_bits = lfsr_sequence(ALTERNATE_TAPS[order], seed, N)
    seq2_bp = to_bipolar(seq2_bits)
    lags_x, xcorr = cross_correlation(seq_bp, seq2_bp)

    ax = axes[1, 0]
    ax.stem(lags_x, xcorr, linefmt="tab:orange", markerfmt="o",
            basefmt="k-")
    ax.set_title("Cross-Correlation (two different 4th-order sequences)")
    ax.set_xlabel("Circular shift (chips)")
    ax.set_ylabel("Normalized correlation")
    ax.grid(True, alpha=0.3)
    ax.annotate("Low & flat: sequences are\nnearly orthogonal",
                xy=(0.45, 0.80), xycoords="axes fraction", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow",
                          alpha=0.8))

    # -- Panel 4: Sequence-length comparison (3rd, 4th, 5th order) --------
    orders = [3, 4, 5]
    colors = ["tab:blue", "tab:green", "tab:red"]
    ax = axes[1, 1]

    for idx, o in enumerate(orders):
        n = 2 ** o - 1
        s = lfsr_sequence(MAXIMAL_LENGTH_TAPS[o], [1] + [0] * (o - 1), n)
        bp = to_bipolar(s)
        ax.step(np.arange(n), bp + idx * 3, where="mid", color=colors[idx],
                linewidth=1.2, label=f"Order {o}:  period = {n}")

    ax.set_title("Sequence Length Grows Exponentially with Order")
    ax.set_xlabel("Chip index")
    ax.set_ylabel("Value (offset for clarity)")
    ax.set_yticks([])
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.3, axis="x")

    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Interactive widget (for Jupyter notebooks)
# ---------------------------------------------------------------------------

def interactive_lfsr():
    """Launch ipywidgets sliders to explore LFSR sequences interactively.

    Sliders:
      - Polynomial order (3-7): sets register size.
      - Tap A / Tap B: bit positions to XOR for the feedback bit.
    The widget shows the bipolar output and its autocorrelation side by side.
    """
    try:
        from ipywidgets import interact, IntSlider, fixed
    except ImportError:
        print("ipywidgets not available -- run in Jupyter for interactive mode")
        return

    def _update(order=4, tap_a=0, tap_b=1):
        # Clamp taps to valid range
        tap_a = min(tap_a, order - 1)
        tap_b = min(tap_b, order - 1)
        taps = sorted(set([tap_a, tap_b]))

        N = 2 ** order - 1
        seed = [1] + [0] * (order - 1)

        try:
            seq = lfsr_sequence(taps, seed, N)
        except Exception:
            print("Invalid tap configuration.")
            return

        bp = to_bipolar(seq)
        lags, acorr = autocorrelation(bp)

        # Check if maximal-length
        unique_states = len(set(tuple(seq[i:i+order])
                                for i in range(len(seq) - order + 1)))
        is_maximal = unique_states == N

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4))

        poly_str = "1"
        for t in sorted(taps):
            exp = t + 1
            if exp == 1:
                poly_str += " + x"
            else:
                poly_str += f" + x^{exp}"
        poly_str += f" + x^{order}"

        ax1.stem(np.arange(N), bp, linefmt="tab:blue", markerfmt=".",
                 basefmt="k-")
        ax1.set_title(f"LFSR Output  |  {poly_str}  |  "
                      f"{'MAXIMAL' if is_maximal else 'non-maximal'}")
        ax1.set_xlabel("Chip index")
        ax1.set_ylabel("Bipolar value")
        ax1.set_ylim(-1.5, 1.5)
        ax1.grid(True, alpha=0.3)

        ax2.stem(lags, acorr, linefmt="tab:green", markerfmt=".",
                 basefmt="k-")
        ax2.set_title("Autocorrelation")
        ax2.set_xlabel("Circular shift")
        ax2.set_ylabel("Normalized correlation")
        ax2.axhline(y=-1/N, color="tab:red", linestyle="--", linewidth=0.8,
                    label=f"-1/N = {-1/N:.3f}")
        ax2.legend(fontsize=8)
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

    interact(
        _update,
        order=IntSlider(value=4, min=3, max=7, step=1,
                        description="Order"),
        tap_a=IntSlider(value=0, min=0, max=6, step=1,
                        description="Tap A"),
        tap_b=IntSlider(value=1, min=0, max=6, step=1,
                        description="Tap B"),
    )


# ---------------------------------------------------------------------------
# Standalone execution
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    matplotlib.use("Agg")
    fig = plot_lfsr_demo()
    fig.savefig("lfsr_demo.png", dpi=150)
    print("Saved lfsr_demo.png")
    plt.show()
