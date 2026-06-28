"""
Phased Arrays — Angle of Arrival & Beam Steering

A phased array is a group of antennas (called "elements") arranged in a
pattern, usually a straight line. Instead of one big antenna, you use many
small ones working together. Why? Because by combining their signals
cleverly, you can electronically "point" the array in any direction — no
moving parts required.

How it works:
  When a radio wave arrives at an angle, it hits each antenna at a slightly
  different time. That time difference shows up as a phase difference (the
  wave is at a different point in its cycle at each element). By measuring
  these phase differences, you can figure out the signal's Angle of Arrival
  (AOA). Conversely, by artificially adding phase shifts before transmitting,
  you can steer the beam in a chosen direction.

Key concepts:
  - Element spacing: typically half a wavelength, which avoids ambiguities
  - Array factor: a mathematical pattern describing how sensitive the array
    is in each direction — like a "sensitivity map"
  - Beam steering: adding a progressive phase shift across elements to move
    the main lobe (peak sensitivity) to a desired angle
  - More elements = narrower main beam = better angular resolution, but the
    array gets physically larger

The math is rooted in the phase delay formula:
  phase_delay[n] = 2 * pi * n * d * sin(theta) / wavelength
where n is the element index, d is the spacing, and theta is the signal
angle measured from boresight (straight ahead).
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SPEED_OF_LIGHT = 3e8  # m/s


# ---------------------------------------------------------------------------
# Core functions (importable by notebook)
# ---------------------------------------------------------------------------

def compute_phase_delays(num_elements, spacing, frequency, angle):
    """Compute the phase delay at each element for a signal arriving at a given angle.

    When a plane wave arrives at an angle, each antenna element receives it at
    a slightly different time. This function calculates the resulting phase
    offset at each element relative to the first element (index 0).

    Parameters
    ----------
    num_elements : int – number of antenna elements in the linear array
    spacing : float – distance between adjacent elements (meters)
    frequency : float – signal frequency (Hz)
    angle : float – angle of arrival in degrees, measured from boresight
                     (0 = straight ahead, positive = right)

    Returns
    -------
    phase_delays : ndarray – phase delay at each element (radians), shape (num_elements,)
    """
    wavelength = SPEED_OF_LIGHT / frequency
    indices = np.arange(num_elements)
    angle_rad = np.deg2rad(angle)
    phase_delays = 2 * np.pi * indices * spacing * np.sin(angle_rad) / wavelength
    return phase_delays


def array_factor(num_elements, spacing, frequency, angles):
    """Compute the array factor (beam pattern) across a range of angles.

    The array factor tells you how sensitive the array is in each direction.
    It is the magnitude of the sum of all element contributions, where each
    element has a phase shift determined by its position and the look angle.
    Think of it as the array's "sensitivity map."

    The result is normalized so the peak equals 1.0.

    Parameters
    ----------
    num_elements : int – number of antenna elements
    spacing : float – element spacing (meters)
    frequency : float – signal frequency (Hz)
    angles : ndarray – angles to evaluate (degrees)

    Returns
    -------
    af : ndarray – normalized array factor magnitude, same shape as angles
    """
    wavelength = SPEED_OF_LIGHT / frequency
    angles_rad = np.deg2rad(angles)
    indices = np.arange(num_elements)

    # Each row is one angle, each column is one element
    # phase_matrix[i, n] = phase contribution of element n at angle i
    phase_matrix = 2 * np.pi * indices[np.newaxis, :] * spacing * \
        np.sin(angles_rad[:, np.newaxis]) / wavelength

    # Sum across elements, take magnitude
    af = np.abs(np.sum(np.exp(1j * phase_matrix), axis=1))

    # Normalize to peak = 1
    af_max = np.max(af)
    if af_max > 0:
        af /= af_max

    return af


def estimate_aoa(phase_delays, spacing, frequency):
    """Estimate the Angle of Arrival from measured phase delays.

    This is the inverse of compute_phase_delays. Given the phase differences
    between adjacent elements, it recovers the signal angle using:
      theta = arcsin(delta_phase * wavelength / (2 * pi * d))

    For a clean signal with half-wavelength spacing, this gives an
    unambiguous answer between -90 and +90 degrees.

    Parameters
    ----------
    phase_delays : ndarray – phase delay at each element (radians)
    spacing : float – element spacing (meters)
    frequency : float – signal frequency (Hz)

    Returns
    -------
    aoa_deg : float – estimated angle of arrival (degrees)
    """
    wavelength = SPEED_OF_LIGHT / frequency

    # Use the average phase difference between consecutive elements for a
    # more robust estimate (reduces noise compared to using just one pair)
    consecutive_diffs = np.diff(phase_delays)
    avg_phase_diff = np.mean(consecutive_diffs)

    # Solve for angle: sin(theta) = avg_phase_diff * wavelength / (2*pi*d)
    sin_theta = avg_phase_diff * wavelength / (2 * np.pi * spacing)

    # Clip to valid arcsin range to handle numerical noise
    sin_theta = np.clip(sin_theta, -1.0, 1.0)
    aoa_deg = np.rad2deg(np.arcsin(sin_theta))

    return float(aoa_deg)


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def plot_array_demo(frequency=1e9, num_elements=8, signal_angle=30):
    """Generate the four-panel overview figure.

    Panels:
      1. Array geometry with incoming wavefront arrows
      2. Array factor beam pattern (polar plot)
      3. Phase delays across elements for AOA estimation
      4. Beam pattern comparison for different array sizes
    """
    wavelength = SPEED_OF_LIGHT / frequency
    spacing = wavelength / 2  # half-wavelength spacing

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Phased Array Fundamentals", fontsize=14, fontweight="bold")

    # --- Panel 1: Array geometry with wavefront arrows ---
    ax = axes[0, 0]
    element_positions = np.arange(num_elements) * spacing

    # Center the array at zero for cleaner visualization
    center_offset = element_positions[-1] / 2
    element_positions_centered = element_positions - center_offset

    # Plot antenna elements
    ax.plot(element_positions_centered / wavelength, np.zeros(num_elements),
            "s", color="tab:blue", markersize=10, label="Antenna elements")

    # Draw incoming wavefront arrows
    angle_rad = np.deg2rad(signal_angle)
    arrow_len = 0.6
    dx = -arrow_len * np.sin(angle_rad)
    dy = -arrow_len * np.cos(angle_rad)
    for pos in element_positions_centered / wavelength:
        ax.annotate("", xy=(pos, 0.05),
                     xytext=(pos - dx, 0.05 - dy),
                     arrowprops=dict(arrowstyle="->", color="tab:red",
                                     lw=1.5))

    # Draw a wavefront line (perpendicular to the arrows)
    wf_half = (num_elements * spacing / wavelength) * 0.6
    wf_center_x = -dx * 0.5
    wf_center_y = 0.05 - dy * 0.5
    perp_dx = np.cos(angle_rad) * wf_half
    perp_dy = -np.sin(angle_rad) * wf_half
    ax.plot([wf_center_x - perp_dx, wf_center_x + perp_dx],
            [wf_center_y - perp_dy, wf_center_y + perp_dy],
            "--", color="tab:red", alpha=0.6, lw=1.5, label="Wavefront")

    ax.set_title(f"Array Geometry (signal at {signal_angle}°)")
    ax.set_xlabel("Position (λ)")
    ax.set_ylabel("")
    ax.set_ylim(-0.5, 1.5)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, alpha=0.3)
    ax.set_aspect("equal")

    # --- Panel 2: Array factor polar plot ---
    ax_polar = fig.add_subplot(2, 2, 2, projection="polar")
    axes[0, 1].remove()  # Remove the rectangular axis, replace with polar

    scan_angles = np.linspace(-90, 90, 721)
    af = array_factor(num_elements, spacing, frequency, scan_angles)
    af_db = 20 * np.log10(af + 1e-12)
    af_db = np.clip(af_db, -40, 0)  # Floor at -40 dB

    # Polar plot: 0 degrees = up (boresight), convert to polar convention
    theta_polar = np.deg2rad(scan_angles)
    ax_polar.plot(theta_polar, af_db + 40, color="tab:blue", lw=1.5)
    ax_polar.fill_between(theta_polar, 0, af_db + 40, alpha=0.15,
                          color="tab:blue")
    ax_polar.set_theta_zero_location("N")
    ax_polar.set_theta_direction(-1)
    ax_polar.set_thetamin(-90)
    ax_polar.set_thetamax(90)
    ax_polar.set_title(f"Array Factor ({num_elements} elements)", pad=20)
    ax_polar.set_rticks([0, 10, 20, 30, 40])
    ax_polar.set_yticklabels(["-40", "-30", "-20", "-10", "0 dB"], fontsize=7)

    # --- Panel 3: Phase delays across elements (AOA estimation) ---
    ax = axes[1, 0]
    phase_delays = compute_phase_delays(num_elements, spacing, frequency,
                                        signal_angle)
    estimated_angle = estimate_aoa(phase_delays, spacing, frequency)

    ax.plot(range(num_elements), np.rad2deg(phase_delays), "o-",
            color="tab:green", lw=2, markersize=8, label="Phase delay")

    # Show the linear fit that AOA estimation uses
    avg_slope = np.mean(np.diff(phase_delays))
    fit_line = np.arange(num_elements) * avg_slope
    ax.plot(range(num_elements), np.rad2deg(fit_line), "--",
            color="tab:orange", lw=1.5, alpha=0.8, label="Linear fit")

    ax.set_title(f"Phase Delays → AOA Estimate: {estimated_angle:.1f}°"
                 f" (true: {signal_angle}°)")
    ax.set_xlabel("Element index")
    ax.set_ylabel("Phase delay (degrees)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # --- Panel 4: Beam pattern comparison for different array sizes ---
    ax = axes[1, 1]
    scan_angles_lin = np.linspace(-90, 90, 721)
    for n in [4, 8, 16]:
        af_n = array_factor(n, spacing, frequency, scan_angles_lin)
        af_n_db = 20 * np.log10(af_n + 1e-12)
        af_n_db = np.clip(af_n_db, -40, 0)
        ax.plot(scan_angles_lin, af_n_db, lw=1.5, label=f"{n} elements")

    ax.set_title("Beam Pattern vs Array Size")
    ax.set_xlabel("Angle (degrees)")
    ax.set_ylabel("Normalized gain (dB)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-90, 90)
    ax.set_ylim(-40, 2)

    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Interactive widget (for Jupyter notebook)
# ---------------------------------------------------------------------------

def interactive_array():
    """Launch ipywidgets sliders for exploring phased array behavior.

    Sliders:
      - Number of elements (2 to 32)
      - Signal angle of arrival (-80 to +80 degrees)
      - Frequency (100 MHz to 10 GHz)
    """
    try:
        from ipywidgets import interact, IntSlider, FloatSlider
    except ImportError:
        print("ipywidgets not available — run in Jupyter for interactive mode")
        return

    def _update(num_elements=8, signal_angle=30, freq_ghz=1.0):
        frequency = freq_ghz * 1e9
        wavelength = SPEED_OF_LIGHT / frequency
        spacing = wavelength / 2

        scan_angles = np.linspace(-90, 90, 721)
        af = array_factor(num_elements, spacing, frequency, scan_angles)
        af_db = 20 * np.log10(af + 1e-12)
        af_db = np.clip(af_db, -40, 0)

        phase_delays = compute_phase_delays(num_elements, spacing, frequency,
                                            signal_angle)
        estimated_angle = estimate_aoa(phase_delays, spacing, frequency)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

        # Beam pattern
        ax1.plot(scan_angles, af_db, color="tab:blue", lw=2)
        ax1.axvline(signal_angle, color="tab:red", ls="--", alpha=0.7,
                    label=f"Signal at {signal_angle}°")
        ax1.fill_between(scan_angles, -40, af_db, alpha=0.1, color="tab:blue")
        ax1.set_title(f"Array Factor ({num_elements} elements, "
                      f"{freq_ghz:.1f} GHz)")
        ax1.set_xlabel("Angle (degrees)")
        ax1.set_ylabel("Gain (dB)")
        ax1.set_xlim(-90, 90)
        ax1.set_ylim(-40, 2)
        ax1.legend(fontsize=8)
        ax1.grid(True, alpha=0.3)

        # Phase delays with AOA estimate
        ax2.plot(range(num_elements), np.rad2deg(phase_delays), "o-",
                 color="tab:green", lw=2, markersize=8)
        ax2.set_title(f"Phase Delays → AOA Estimate: "
                      f"{estimated_angle:.1f}° (true: {signal_angle}°)")
        ax2.set_xlabel("Element index")
        ax2.set_ylabel("Phase delay (degrees)")
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

    interact(
        _update,
        num_elements=IntSlider(value=8, min=2, max=32, step=1,
                               description="Elements"),
        signal_angle=FloatSlider(value=30, min=-80, max=80, step=1,
                                 description="Angle (°)"),
        freq_ghz=FloatSlider(value=1.0, min=0.1, max=10.0, step=0.1,
                             description="Freq (GHz)"),
    )


# ---------------------------------------------------------------------------
# Standalone execution
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    matplotlib.use("Agg")
    fig = plot_array_demo()
    fig.savefig("phased_array.png", dpi=150)
    print("Saved phased_array.png")

    plt.show()
