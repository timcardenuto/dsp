"""
GPS Gold Codes & CDMA (Code Division Multiple Access)

GPS satellites all broadcast on the SAME frequency at the SAME time. How does
your phone tell them apart? Each satellite is assigned a unique pseudo-random
bit sequence called a "Gold code" (a type of PRN — Pseudo-Random Noise). These
codes have a special mathematical property: they are nearly uncorrelated with
each other. That means if you multiply a received signal by one satellite's
code, that satellite's data pops out clearly while every other satellite's
signal smears into low-level noise.

How it works step by step:
  1. TRANSMIT — Each satellite multiplies its slow data bits by its fast PRN
     code (1023 chips). This "spreads" the signal across a wide bandwidth.
  2. COMBINE — All satellite signals add together in the air (plus noise).
     The result looks like random noise to anyone without the codes.
  3. RECEIVE — Your receiver multiplies the combined signal by the PRN code
     of the satellite it wants. This "despreads" that satellite's data back
     to its original form. All other satellites remain spread (noise-like)
     and are removed by a simple averaging (low-pass) filter.

This is CDMA in action: multiple users share the same channel simultaneously,
separated only by their unique codes — no time slots, no frequency bands.

Gold code generation uses two 10-stage Linear Feedback Shift Registers (LFSRs)
called G1 and G2. G1 is the same for every satellite; G2's output taps differ
per satellite, producing 1023-chip codes with excellent cross-correlation.
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# GPS PRN tap assignments (SV number -> G2 output taps)
# ---------------------------------------------------------------------------

# Each GPS satellite (Space Vehicle) is identified by a PRN number. The G2
# output taps determine which two register outputs are XORed to produce that
# satellite's unique code. These assignments are defined in the GPS ICD.
GPS_SV_TAPS = {
    1: [2, 6],   2: [3, 7],   3: [4, 8],   4: [5, 9],   5: [1, 9],
    6: [2, 10],  7: [1, 8],   8: [2, 9],   9: [3, 10],  10: [2, 3],
    11: [3, 4],  12: [5, 6],  13: [6, 7],  14: [7, 8],  15: [8, 9],
    16: [9, 10], 17: [1, 4],  18: [2, 5],  19: [3, 6],  20: [4, 7],
    21: [5, 8],  22: [6, 9],  23: [1, 3],  24: [4, 6],  25: [5, 7],
    26: [6, 8],  27: [7, 9],  28: [8, 10], 29: [1, 6],  30: [2, 7],
    31: [3, 8],  32: [4, 9],
}


# ---------------------------------------------------------------------------
# Core functions (importable by notebook)
# ---------------------------------------------------------------------------

def generate_generic_prn(order, seed, feedback_taps, output_taps):
    """Generate a pseudo-random bit sequence using a generic LFSR.

    A Linear Feedback Shift Register (LFSR) is a chain of flip-flops that
    shifts bits on every clock tick. Some bits are XORed together and fed
    back to the input — the choice of which bits (the "taps") determines
    the sequence produced. With the right taps you get a maximal-length
    sequence: 2^order - 1 bits long before it repeats.

    Parameters
    ----------
    order : int
        Number of shift register stages (bits). A 10-stage LFSR produces
        a sequence of length 2^10 - 1 = 1023.
    seed : array-like of int
        Initial fill of the register, length must equal *order*.
        Typically all ones for GPS.
    feedback_taps : list of int
        1-based positions XORed together to create the feedback bit that
        enters the first register on each clock tick.
    output_taps : list of int
        1-based positions XORed together to produce each output bit.
        Use [order] to simply output the last register stage.

    Returns
    -------
    ndarray of int (0s and 1s), length 2^order - 1
    """
    N = 2**order - 1
    X = np.array(seed, dtype=int).copy()
    output = np.zeros(N, dtype=int)

    for i in range(N):
        # Compute feedback bit from the designated taps
        feedback = 0
        for tap in feedback_taps:
            feedback ^= X[tap - 1]         # taps are 1-based

        # Shift register right: move each bit to the next position
        # (highest index = oldest bit)
        X[1:] = X[:-1]
        X[0] = feedback

        # Compute output bit by XORing the output taps
        out_bit = 0
        for tap in output_taps:
            out_bit ^= X[tap - 1]
        output[i] = out_bit

    return output


def generate_gps_prn(sv_number):
    """Generate the 1023-chip Gold code for a GPS satellite.

    Each GPS satellite broadcasts a unique Gold code built from two LFSRs:
      - G1: feedback polynomial x^10 + x^3 + 1  (taps 3 and 10)
      - G2: feedback polynomial x^10 + x^9 + x^8 + x^6 + x^3 + x^2 + 1
    Both start with all-ones. The satellite-specific part is which two G2
    register stages are XORed to form the G2 output. The final PRN chip is
    G1_output XOR G2_output on each clock cycle.

    Parameters
    ----------
    sv_number : int
        Satellite PRN number (1-32).

    Returns
    -------
    ndarray of int (0s and 1s), length 1023
    """
    if sv_number not in GPS_SV_TAPS:
        raise ValueError(f"SV number must be 1-32, got {sv_number}")

    order = 10
    seed = np.ones(order, dtype=int)

    # G1 LFSR — same for every satellite
    # Feedback taps: registers 3 and 10; output: register 10
    g1 = generate_generic_prn(order, seed, [3, 10], [10])

    # G2 LFSR — same feedback polynomial, but output taps vary per satellite
    g2_taps = GPS_SV_TAPS[sv_number]
    g2 = generate_generic_prn(order, seed, [2, 3, 6, 8, 9, 10], g2_taps)

    # Gold code = G1 XOR G2
    prn = np.bitwise_xor(g1, g2)
    return prn


def gps_transmit(messages, prn_numbers):
    """Simulate GPS transmission from multiple satellites.

    Each satellite's message is "spread" by multiplying it sample-by-sample
    with a bipolar (+1/-1) version of its PRN code. All spread signals are
    then summed together — just as they combine in the air.

    Parameters
    ----------
    messages : list of ndarray
        Baseband message signal for each satellite (same length).
    prn_numbers : list of int
        PRN (SV) number for each satellite.

    Returns
    -------
    combined : ndarray
        Sum of all spread signals (what a receiver antenna picks up).
    spread_signals : list of ndarray
        Individual spread signal for each satellite.
    bipolar_prns : list of ndarray
        Bipolar PRN codes repeated/trimmed to message length.
    """
    n_samples = len(messages[0])
    spread_signals = []
    bipolar_prns = []

    for msg, sv in zip(messages, prn_numbers):
        prn = generate_gps_prn(sv)

        # Map 0/1 to -1/+1 (bipolar) — required so that multiplying twice
        # gives back +1 (despreading property: code * code = 1)
        bipolar = 2 * prn.astype(float) - 1

        # Tile the 1023-chip code to cover the entire message length
        repeats = int(np.ceil(n_samples / len(bipolar)))
        bipolar_long = np.tile(bipolar, repeats)[:n_samples]

        spread = msg * bipolar_long
        spread_signals.append(spread)
        bipolar_prns.append(bipolar_long)

    combined = np.sum(spread_signals, axis=0)
    return combined, spread_signals, bipolar_prns


def gps_receive(signal, prn_number, averaging_length=31):
    """Simulate GPS reception — extract one satellite's data from the mix.

    Two steps:
      1. DESPREAD — multiply the received signal by the target satellite's
         PRN code. Because code * code = +1, the target data is restored.
         Other satellites' data remains spread (noise-like).
      2. FILTER — a simple moving-average (low-pass) filter smooths out the
         remaining spread-spectrum noise, revealing the original message.

    Parameters
    ----------
    signal : ndarray
        Received combined signal (output of gps_transmit).
    prn_number : int
        Which satellite to extract (1-32).
    averaging_length : int
        Number of samples for the moving-average filter. Larger values
        give cleaner output but introduce more delay.

    Returns
    -------
    despread : ndarray
        Signal after multiplying by the PRN code (before filtering).
    filtered : ndarray
        Signal after low-pass averaging (recovered message).
    """
    prn = generate_gps_prn(prn_number)
    bipolar = 2 * prn.astype(float) - 1
    n_samples = len(signal)
    repeats = int(np.ceil(n_samples / len(bipolar)))
    bipolar_long = np.tile(bipolar, repeats)[:n_samples]

    # Step 1: despread
    despread = signal * bipolar_long

    # Step 2: moving-average low-pass filter
    kernel = np.ones(averaging_length) / averaging_length
    filtered = np.convolve(despread, kernel, mode='same')

    return despread, filtered


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def plot_gps_demo():
    """Generate the GPS TX/RX demo across two figures (6+ panels).

    Figure 1 — Transmit side:
      Shows how three satellites' signals are created and combined.
    Figure 2 — Receive side:
      Shows how the receiver pulls one satellite's message back out.
    """
    # -- Build three satellite messages (simple sine waves) ----------------
    fs = 10                                  # symbol rate
    T = 1 / fs
    chips_per_symbol = 31                    # PRN chips per message sample
    n_symbols = 100
    n_samples = n_symbols * chips_per_symbol
    t = np.linspace(0, n_symbols * T, n_samples, endpoint=False)

    msg1 = np.sin(2 * np.pi * fs * t)       # SV 1 message
    msg2 = np.sin(2 * np.pi * fs * 0.7 * t) # SV 5 message (different freq)
    msg3 = np.sin(2 * np.pi * fs * 1.3 * t) # SV 16 message

    messages = [msg1, msg2, msg3]
    prn_numbers = [1, 5, 16]

    combined, spread_signals, bipolar_prns = gps_transmit(
        messages, prn_numbers
    )

    # -- Figure 1: Transmit -----------------------------------------------
    fig1, axes1 = plt.subplots(4, 1, figsize=(12, 10))
    fig1.suptitle("GPS CDMA — Transmit Side", fontsize=14)

    ax = axes1[0]
    ax.plot(t, msg1, color="tab:blue")
    ax.set_title("Original Message (Satellite PRN 1)")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Amplitude")
    ax.grid(True, alpha=0.3)

    ax = axes1[1]
    ax.plot(t, bipolar_prns[0], color="tab:orange", linewidth=0.5)
    ax.set_title("PRN 1 Gold Code (1023-chip spreading sequence, bipolar)")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Amplitude")
    ax.set_ylim(-1.5, 1.5)
    ax.grid(True, alpha=0.3)

    ax = axes1[2]
    ax.plot(t, spread_signals[0], color="tab:green", linewidth=0.5)
    ax.set_title("Spread Signal = Message x PRN Code (looks like noise)")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Amplitude")
    ax.grid(True, alpha=0.3)

    ax = axes1[3]
    ax.plot(t, combined, color="tab:purple", linewidth=0.5)
    ax.set_title(
        f"Combined RF Signal (PRNs {prn_numbers} summed — "
        "what the antenna receives)"
    )
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Amplitude")
    ax.grid(True, alpha=0.3)

    fig1.tight_layout()

    # -- Figure 2: Receive -------------------------------------------------
    despread, filtered = gps_receive(combined, prn_number=1,
                                     averaging_length=chips_per_symbol)

    fig2, axes2 = plt.subplots(3, 1, figsize=(12, 8))
    fig2.suptitle("GPS CDMA — Receive Side (extracting PRN 1)", fontsize=14)

    ax = axes2[0]
    ax.plot(t, combined, color="tab:purple", linewidth=0.5)
    ax.set_title("Received Combined Signal (all satellites mixed together)")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Amplitude")
    ax.grid(True, alpha=0.3)

    ax = axes2[1]
    ax.plot(t, despread, color="tab:red", linewidth=0.5)
    ax.set_title("After Despreading with PRN 1 (multiply by same code)")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Amplitude")
    ax.grid(True, alpha=0.3)

    ax = axes2[2]
    ax.plot(t, filtered, color="tab:red", linewidth=1.5,
            label="Recovered message")
    ax.plot(t, msg1, color="tab:blue", linewidth=1.5, alpha=0.6,
            linestyle="--", label="Original message")
    ax.set_title("Filtered/Recovered Message vs Original")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Amplitude")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    fig2.tight_layout()
    return fig1, fig2


def plot_correlation_demo(sv_list=None):
    """Show cross-correlation between GPS Gold codes.

    Gold codes are designed so that any two different codes have very low
    cross-correlation. This plot shows a matrix: each cell is the normalized
    peak cross-correlation between two PRN codes. The diagonal (a code
    correlated with itself) is 1.0; off-diagonal values should be very
    small — proving the codes are nearly orthogonal.

    Parameters
    ----------
    sv_list : list of int or None
        Which PRN numbers to compare. Defaults to [1..8].

    Returns
    -------
    fig : matplotlib Figure
    """
    if sv_list is None:
        sv_list = list(range(1, 9))

    n = len(sv_list)
    codes = []
    for sv in sv_list:
        prn = generate_gps_prn(sv)
        bipolar = 2 * prn.astype(float) - 1
        codes.append(bipolar)

    # Build correlation matrix
    corr_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            # Normalized cross-correlation at zero lag
            corr_matrix[i, j] = (
                np.dot(codes[i], codes[j])
                / np.sqrt(np.dot(codes[i], codes[i]) * np.dot(codes[j], codes[j]))
            )

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle("GPS Gold Code Cross-Correlation", fontsize=14)

    # Heatmap
    ax = axes[0]
    im = ax.imshow(corr_matrix, cmap="RdBu_r", vmin=-0.1, vmax=1.0,
                   interpolation="nearest")
    ax.set_xticks(range(n))
    ax.set_xticklabels([f"PRN {s}" for s in sv_list], rotation=45, ha="right")
    ax.set_yticks(range(n))
    ax.set_yticklabels([f"PRN {s}" for s in sv_list])
    ax.set_title("Correlation Matrix (zero-lag)")
    fig.colorbar(im, ax=ax, shrink=0.8)

    # Annotate cells with values
    for i in range(n):
        for j in range(n):
            color = "white" if corr_matrix[i, j] > 0.5 else "black"
            ax.text(j, i, f"{corr_matrix[i, j]:.3f}", ha="center",
                    va="center", fontsize=7, color=color)

    # Full cross-correlation example between two codes
    ax = axes[1]
    c1, c2 = codes[0], codes[1]
    full_corr = np.correlate(c1, c2, mode='full') / len(c1)
    lags = np.arange(-len(c1) + 1, len(c1))
    ax.plot(lags, full_corr, color="tab:blue", linewidth=0.5)
    ax.set_title(
        f"Full Cross-Correlation: PRN {sv_list[0]} vs PRN {sv_list[1]}"
    )
    ax.set_xlabel("Lag (chips)")
    ax.set_ylabel("Normalized Correlation")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Interactive widget (for Jupyter notebook)
# ---------------------------------------------------------------------------

def interactive_gps():
    """Launch ipywidgets sliders to explore GPS CDMA interactively.

    Controls:
      - Number of satellites transmitting (1-6)
      - Which PRN number the receiver locks onto
    """
    try:
        from ipywidgets import interact, IntSlider
    except ImportError:
        print("ipywidgets not available — run in Jupyter for interactive mode")
        return

    available_svs = [1, 5, 7, 10, 16, 23]

    def _update(n_satellites=3, receive_prn=1):
        svs_in_use = available_svs[:n_satellites]
        if receive_prn not in svs_in_use:
            receive_prn = svs_in_use[0]

        # Build messages
        fs = 10
        T = 1 / fs
        chips_per_symbol = 31
        n_symbols = 100
        n_samples = n_symbols * chips_per_symbol
        t = np.linspace(0, n_symbols * T, n_samples, endpoint=False)

        messages = []
        for idx, sv in enumerate(svs_in_use):
            freq_scale = 1.0 + 0.3 * idx
            messages.append(np.sin(2 * np.pi * fs * freq_scale * t))

        combined, _, _ = gps_transmit(messages, svs_in_use)
        despread, filtered = gps_receive(combined, receive_prn,
                                         averaging_length=chips_per_symbol)

        # Find the original message for the selected PRN
        rx_idx = svs_in_use.index(receive_prn)
        original = messages[rx_idx]

        fig, axes = plt.subplots(3, 1, figsize=(12, 8))
        fig.suptitle(
            f"GPS CDMA — {n_satellites} satellites, receiving PRN {receive_prn}",
            fontsize=14,
        )

        ax = axes[0]
        ax.plot(t, combined, color="tab:purple", linewidth=0.5)
        ax.set_title(
            f"Combined Signal (SVs {svs_in_use})"
        )
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Amplitude")
        ax.grid(True, alpha=0.3)

        ax = axes[1]
        ax.plot(t, despread, color="tab:red", linewidth=0.5)
        ax.set_title(f"After Despreading with PRN {receive_prn}")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Amplitude")
        ax.grid(True, alpha=0.3)

        ax = axes[2]
        ax.plot(t, filtered, color="tab:red", linewidth=1.5,
                label="Recovered")
        ax.plot(t, original, color="tab:blue", linewidth=1.5, alpha=0.6,
                linestyle="--", label="Original")
        ax.set_title("Recovered vs Original Message")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Amplitude")
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

    interact(
        _update,
        n_satellites=IntSlider(value=3, min=1, max=6, step=1,
                               description="Satellites"),
        receive_prn=IntSlider(value=1, min=1, max=23, step=1,
                              description="Receive PRN"),
    )


# ---------------------------------------------------------------------------
# Standalone execution
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    matplotlib.use("Agg")

    fig1, fig2 = plot_gps_demo()
    fig1.savefig("gps_tx_demo.png", dpi=150)
    print("Saved gps_tx_demo.png")
    fig2.savefig("gps_rx_demo.png", dpi=150)
    print("Saved gps_rx_demo.png")

    fig3 = plot_correlation_demo()
    fig3.savefig("gps_correlation_demo.png", dpi=150)
    print("Saved gps_correlation_demo.png")

    plt.show()
