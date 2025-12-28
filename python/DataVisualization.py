from pickle import NONE
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from scipy import signal
from scipy.optimize import curve_fit

# Teledyne Lecroy WaveAce 2024: 200 MHz, 2 GS/s
SAMPLE_RATE_GS = 2.0  # GigaSamples per second
NS_PER_SAMPLE = 1.0 / SAMPLE_RATE_GS  # 0.5 ns per sample
SAMPLES_PER_NS = 1.0 / NS_PER_SAMPLE  # 2.0 samples per ns

def ns_to_samples(ns):
    """Convert nanoseconds to sample indices"""
    return int(ns * SAMPLES_PER_NS)

# then the 350 samples we aare interesed in correspond to 175 ns
def samples_to_ns(samples):
    """Convert sample indices to nanoseconds"""
    return samples / SAMPLES_PER_NS

with open('5_mit_waveform_data.txt', 'r') as f:
    lines = f.readlines()
    waveforms = []

    for i in range(0, len(lines), 3):
        line = lines[i].strip()
        if not line:
            continue

        if line.startswith("CH1:"):
            line = line[4:]

        data = [x for x in line.split(',') if x]
        data = data[:350]
        waveform = [(float(x) - 255 if float(x) > 128 else float(x)) for x in data]
        waveforms.append(waveform)

print(f"Loaded {len(waveforms)} waveforms from file")

# Do signal preprocessing-...
waveforms_clean = []
for waveform_orig in waveforms:
    waveform_inverted = [-x for x in waveform_orig]
    waveforms_clean.append(waveform_inverted)

# Use first 50 samples for baseline estimation (25 ns @ 2 GS/s)
waveforms_baseline_corrected = []
baselines = []
max_amplitudes = []
rise_times = []
max_indices = []

for waveform in waveforms_clean:
    # Baseline correction using first 50 samples
    baseline = np.mean(waveform[:50])
    baselines.append(baseline)
    #waveform_baseline_corrected = [x - baseline for x in waveform]
    waveform_baseline_corrected = [max(0, x - baseline) for x in waveform] # were clipping to prevent negative values
    waveforms_baseline_corrected.append(waveform_baseline_corrected)

    # Find peak amplitude and index
    max_amp = max(waveform_baseline_corrected)
    max_amplitudes.append(max_amp)
    max_idx = np.argmax(waveform_baseline_corrected)
    max_indices.append(max_idx)

    # Calculate rise time, we'll try 10% to 90% of peak
    threshold_10 = 0.1 * max_amp
    threshold_90 = 0.9 * max_amp
    idx_10 = next((i for i, v in enumerate(waveform_baseline_corrected) if v > threshold_10), None)
    idx_90 = next((i for i, v in enumerate(waveform_baseline_corrected) if v > threshold_90), None)
    #if idx_10 and idx_90:
    if idx_10 is not None and idx_90 is not None:
        rise_times.append(idx_90 - idx_10)
    else:
        rise_times.append(0)


def detect_pileup(waveform, max_idx, qtot, psd_factor, threshold_factor=0.25):
    """
    Detect pile-up events using PSD plot analysis.
    Reference: https://uu.diva-portal.org/smash/get/diva2:925705/FULLTEXT01.pdf Site 45-48
    Binda's method: Events above the neutron cluster in the PSD plot are rejected as pile-up.
    This is more reliable than peak detection for NE213 since pile-up distorts both
    amplitude and decay characteristics.

    Key indicators:
    - Abnormally high PSD factor for given energy
    - Events that fall outside expected neutron/gamma clusters

    Args:
        waveform: Baseline-corrected waveform
        max_idx: Index of maximum amplitude
        qtot: Total integrated charge
        psd_factor: Calculated PSD parameter (Qtot-Qshort)/Qtot
        threshold_factor: Relative threshold for outlier detection

    Returns:
        True if pile-up detected, False otherwise
    """
    if max_idx < 20 or max_idx >= len(waveform) - 10:
        return True

    # Check for abnormally high PSD factor (above typical neutron range)
    # For NE213, neutron PSD typically < 0.35, pile-up can push this higher
    if psd_factor > 0.5:  #TODO AdrGos: Verify with Felix if this threshold is valid for NE213 or better to use 0.4-0.45 to reject more pile-up
        return True

    # Checks for double peaks in rising edge (20 samples before peak)
    rising_edge = waveform[max(0, max_idx-20):max_idx]
    if len(rising_edge) > 10:
        peaks, _ = signal.find_peaks(rising_edge, height=threshold_factor*max(rising_edge))
        if len(peaks) > 1:
            return True

    # Checks for abnormal decay e.g sudden increases after pea
    if max_idx + 30 < len(waveform):
        decay_region = waveform[max_idx:max_idx+30]
        # Check if any point increases by more than 15% after peak
        for i in range(1, len(decay_region)-1):
            if decay_region[i] > decay_region[i-1] * 1.15:
                return True

    return False


def calculate_psd_parameters(waveforms, max_indices, short_gate_ns=25, total_gate_ns=150,
                             pre_trigger_ns=10):
    """
    Calculate PSD parameters using Binda's Charge Comparison Method (Equation 8).

    Binda Equation (8): PSD = (Qtot - Qshort) / Qtot
    Reference: https://uu.diva-portal.org/smash/get/diva2:925705/FULLTEXT01.pdf Site 31/32 Figgure 8

    Where:
    - Qtot: Total charge integration (pre_trigger to total_gate)
    - Qshort: Short gate integration (pre_trigger to short_gate)
    - Higher PSD indicates neutrons (longer tail)
    - Lower PSD indicates gammas (shorter tail)

    Args:
        waveforms: List of baseline-corrected waveforms
        max_indices: List of peak position indices
        short_gate_ns: Duration of short integration gate (ns)
        total_gate_ns: Duration of total integration gate (ns)
        pre_trigger_ns: Time before peak to start integration (ns)

    Returns:
        psd_values: PSD parameter (Qtot-Qshort)/Qtot for each pulse
        qtot_values: Total integrated charge
        qshort_values: Short gate charge
        valid_flags: Boolean flags for valid events (pile-up rejected)
    """
    short_gate_samples = ns_to_samples(short_gate_ns)
    total_gate_samples = ns_to_samples(total_gate_ns)
    pre_trigger_samples = ns_to_samples(pre_trigger_ns)

    psd_values = []
    qtot_values = []
    qshort_values = []
    valid_flags = []

    for i, (waveform, max_idx) in enumerate(zip(waveforms, max_indices)):
        if max_idx is None or max_idx >= len(waveform) - 10:
            psd_values.append(0)
            qtot_values.append(0)
            qshort_values.append(0)
            valid_flags.append(False)
            continue

        # Both gates start from the same point (beginning of pulse rise)
        start_idx = max(0, max_idx - pre_trigger_samples)

        # Qtot: Total charge integration (captures full pulse including slow component)
        total_end = min(len(waveform), start_idx + total_gate_samples)
        qtot = np.sum(waveform[start_idx:total_end])

        # Qshort: Short gate integration (captures mainly fast component)
        short_end = min(len(waveform), start_idx + short_gate_samples)
        qshort = np.sum(waveform[start_idx:short_end])

        # PSD parameter calculation form binda
        # The slow component fraction = (Qtot - Qshort) / Qtot
        # Neutrons have longer decay (more slow component) -> higher PSD
        # Gammas have shorter decay (less slow component) -> lower PSD
        if qtot > 0:
            psd = (qtot - qshort) / qtot
        else:
            psd = 0

        # Pile-up detection using PSD-based method (Binda Section 8.1.2)
        is_pileup = detect_pileup(waveform, max_idx, qtot, psd)

        psd_values.append(psd)
        qtot_values.append(qtot)
        qshort_values.append(qshort)
        valid_flags.append(not is_pileup)

    return psd_values, qtot_values, qshort_values, valid_flags


def calculate_fom(psd_values, qtot_values, energy_threshold=10, valid_flags=None):
    """
    Calculate Figure of Merit (FOM) for neutron-gamma separation.

    From Binda Equation (9):
    FOM = DELTA_PEAK / (FWHM_gamma + FWHM_neutron)

    Reference: https://uu.diva-portal.org/smash/get/diva2:925705/FULLTEXT01.pdf Site 34 Table 4
    - FOM > 2.0: Excellent discrimination (Bindas JET results)
    - FOM > 1.0: Good discrimination
    - FOM < 1.0: Poor discrimination

    Note AdrGos: Binda achieved FOM values of 2.25-2.35 for NE213 at JET
    Reference: https://www.mdpi.com/2076-3417/14/13/5532 Site 3 - 12

    FIXED: Better peak finding and FWHM calculation matching Binda's method

    Args:
        psd_values: PSD parameters (Qtot-Qshort)/Qtot
        qtot_values: Total charges (energy proxy)
        energy_threshold: Minimum energy for discrimination
        valid_flags: Boolean flags for valid events

    Returns:
        fom: Figure of merit
        gamma_peak: Peak position for gamma events
        neutron_peak: Peak position for neutron events
    """
    # Filter by energy threshold and validity
    if valid_flags is not None:
        valid_idx = [i for i, (e, v) in enumerate(zip(qtot_values, valid_flags))
                     if e > energy_threshold and v]
    else:
        valid_idx = [i for i, e in enumerate(qtot_values) if e > energy_threshold]

    if len(valid_idx) < 20:
        return 0, 0, 0

    filtered_psd = [psd_values[i] for i in valid_idx]

    # Create histogram for peak finding (Binda Figure 16)
    hist, bin_edges = np.histogram(filtered_psd, bins=100)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Find peaks with proper prominence to avoid noise
    peaks, _ = signal.find_peaks(hist, height=max(hist)*0.15, prominence=max(hist)*0.1, distance=10)

    if len(peaks) < 2:
        return 0, 0, 0

    # Sort peaks by position (lower PSD = gamma, higher PSD = neutron)
    peaks = sorted(peaks)
    gamma_peak_idx = peaks[0]
    neutron_peak_idx = peaks[-1]

    gamma_peak = bin_centers[gamma_peak_idx]
    neutron_peak = bin_centers[neutron_peak_idx]

    separation = abs(neutron_peak - gamma_peak)

    # site  23 figure 7
    def calculate_fwhm(hist, bin_centers, peak_idx):
        peak_height = hist[peak_idx]
        half_max = peak_height / 2.0

        # Find left and right half-maximum points
        left_idx = peak_idx
        while left_idx > 0 and hist[left_idx] > half_max:
            left_idx -= 1

        right_idx = peak_idx
        while right_idx < len(hist)-1 and hist[right_idx] > half_max:
            right_idx += 1

        if left_idx < peak_idx < right_idx:
            fwhm = bin_centers[right_idx] - bin_centers[left_idx]
        else:
            # Fallback to std method fails..
            fwhm = 0.01

        return max(fwhm, 0.001)  # avoid division byy zero

    gamma_fwhm = calculate_fwhm(hist, bin_centers, gamma_peak_idx)
    neutron_fwhm = calculate_fwhm(hist, bin_centers, neutron_peak_idx)

    fom = separation / (gamma_fwhm + neutron_fwhm)

    return fom, gamma_peak, neutron_peak


# Do intial analysis...
# TODO AdrGos: Check these gate values with Felix again
psd, qtot, qshort, valid_flags = calculate_psd_parameters(
    waveforms_baseline_corrected, max_indices,
    short_gate_ns=25,
    total_gate_ns=150,
    pre_trigger_ns=10
)


# Figure 1: Plot Raw Waveforms...
fig1 = plt.figure("Raw Waveforms - NE213 Analysis", figsize=(14, 5))
ax1 = fig1.add_subplot(121)
time_axis = np.arange(350) * NS_PER_SAMPLE  # Time axis in nanoseconds

# Plot all waveforms...
for i, waveform in enumerate(waveforms_baseline_corrected):
    ax1.plot(time_axis[:len(waveform)], waveform, alpha=0.3, linewidth=0.5)

ax1.set_xlabel('Time (ns)')
ax1.set_ylabel('Amplitude (baseline corrected)')
ax1.set_title('NE213 Scintillation Pulses (D-D Fusion)')
ax1.grid(True, alpha=0.3)
ax1.axvline(x=25, color='green', linestyle='--', alpha=0.5, label='Short gate end (25 ns)')
ax1.axvline(x=150, color='red', linestyle='--', alpha=0.5, label='Total gate end (150 ns)')
ax1.legend()

# Plot Average pulse shapes comparison from binda...
ax2 = fig1.add_subplot(122)
sorted_by_psd = sorted(zip(psd, waveforms_baseline_corrected, valid_flags), key=lambda x: x[0])
# Filter out invalid events
valid_sorted = [(psd_val, wf) for psd_val, wf, valid in sorted_by_psd if valid]

if len(valid_sorted) > 10:
    low_psd = [wf for psd_val, wf in valid_sorted[:len(valid_sorted)//3]]
    high_psd = [wf for psd_val, wf in valid_sorted[-len(valid_sorted)//3:]]

    if low_psd and high_psd:
        avg_low = np.mean(low_psd, axis=0)
        avg_high = np.mean(high_psd, axis=0)
        ax2.plot(time_axis[:len(avg_low)], avg_low, 'b-', label='Low PSD (Y-Gamma)', linewidth=2)
        ax2.plot(time_axis[:len(avg_high)], avg_high, 'r-', label='High PSD (n-Neutron)', linewidth=2)
        ax2.set_xlabel('Time (ns)')
        ax2.set_ylabel('Amplitude')
        ax2.set_title('Average Pulse Shapes: Neutron vs Gamma (log scale)')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_yscale('log')
        ax2.set_ylim(bottom=0.1)
        ax2.axvline(x=25, color='green', linestyle='--', alpha=0.3)
        ax2.axvline(x=150, color='red', linestyle='--', alpha=0.3)

plt.tight_layout()

# Figure 2: crreates 2D PSD Plot
fig2, ax_psd = plt.subplots(figsize=(10, 8), num="PSD Analysis: NE213 Neutron-Gamma Discrimination")
# Filter for valid events
valid_qtot = [qtot[i] for i, v in enumerate(valid_flags) if v]
valid_psd = [psd[i] for i, v in enumerate(valid_flags) if v]

scatter = ax_psd.scatter(valid_qtot, valid_psd, c=valid_psd, cmap='coolwarm', alpha=0.6, s=20, edgecolors='black', linewidth=0.5)
ax_psd.set_xlabel('Qtot - Total Charge (Energy Proxy)', fontsize=12)
ax_psd.set_ylabel('PSD = (Qtot - Qshort) / Qtot', fontsize=12)
ax_psd.set_title('NE213 Pulse Shape Discrimination (Binda Method)\n' + 'D-D Fusion\n' + '(Neutrons: high PSD, Gammas: low PSD)', fontsize=14)
ax_psd.grid(True, alpha=0.3)
cbar = plt.colorbar(scatter, ax=ax_psd)
cbar.set_label('PSD Parameter', fontsize=10)

# Add discrimination line
if len(valid_psd) > 0:
    median_psd = np.median(valid_psd)
    ax_psd.axhline(y=median_psd, color='green', linestyle='--', linewidth=2, label=f'Discrimination threshold (median PSD={median_psd:.3f})')
    ax_psd.legend()

# Figure 3: Interactive Histogram with Parameter Control
fig3 = plt.figure("Interactive PSD Histogram - NE213 Optimization", figsize=(12, 8))
plt.subplots_adjust(bottom=0.35)
ax_hist = fig3.add_subplot(111)

# Initial parameters in NANOSECONDS
# TODO AdrGGos: Maybe check these again with felix
SHORT_GATE_NS_INIT = 25
TOTAL_GATE_NS_INIT = 150
ENERGY_THRESHOLD_INIT = 10
BINS_INIT = 100

def update_plot(short_gate_ns, total_gate_ns, energy_threshold, nbins):
    """
    Update PSD histogram with new integration parameters.

    This allows optimization of integration windows for your specific:
    - Detector geometry and volume
    - PMT characteristics
    - Neutron energy range (D-D at 2.5 MeV)
    - Background gamma spectrum
    """
    psd_new, qtot_new, _, valid_new = calculate_psd_parameters(
        waveforms_baseline_corrected, max_indices,
        short_gate_ns=short_gate_ns,
        total_gate_ns=total_gate_ns,
        pre_trigger_ns=10
    )

    # Filter by energy and validity
    valid_idx = [i for i, (e, v) in enumerate(zip(qtot_new, valid_new))
                 if e > energy_threshold and v]
    filtered_psd = [psd_new[i] for i in valid_idx]

    ax_hist.clear()

    if len(filtered_psd) > 0:
        _, _, _ = ax_hist.hist(filtered_psd, bins=int(nbins), color='blue', alpha=0.7, edgecolor='black')
        fom, gamma_pk, neutron_pk = calculate_fom(psd_new, qtot_new, energy_threshold, valid_new)

        ax_hist.set_xlabel('PSD = (Qtot - Qshort) / Qtot', fontsize=12)
        ax_hist.set_ylabel('Counts', fontsize=12)

        short_gate_samples = ns_to_samples(short_gate_ns)
        total_gate_samples = ns_to_samples(total_gate_ns)

        title = f'NE213 PSD Distribution (FOM={fom:.3f})\n'
        title += f'Valid events: {len(filtered_psd)} | Short gate: {short_gate_ns:.1f} ns ({short_gate_samples} samples) | '
        title += f'Total gate: {total_gate_ns:.1f} ns ({total_gate_samples} samples)'

        if fom > 0:
            title += f'\nY peak={gamma_pk:.3f} | n peak={neutron_pk:.3f}'
            if fom > 2.0:
                title += ' | VERY GGOOD (matching Binda JET results)'
            elif fom > 1.0:
                title += ' | GOOD discrimination'
            elif fom > 0.5:
                title += ' | MODERATE discrimination'
            else:
                title += ' | POOR - optimize windows'

        ax_hist.set_title(title, fontsize=11)
        ax_hist.grid(True, alpha=0.3)

    fig3.canvas.draw_idle()

# Create all sliders...
ax_short_gate = plt.axes([0.2, 0.25, 0.65, 0.03])
slider_short_gate = Slider(ax_short_gate, 'Short Gate (ns)', 10, 100, valinit=SHORT_GATE_NS_INIT, valstep=0.5)

ax_total_gate = plt.axes([0.2, 0.20, 0.65, 0.03])
slider_total_gate = Slider(ax_total_gate, 'Total Gate (ns)', 50, 300, valinit=TOTAL_GATE_NS_INIT, valstep=1)

ax_energy = plt.axes([0.2, 0.15, 0.65, 0.03])
slider_energy = Slider(ax_energy, 'Qtot Threshold (energy)', 0, 100, valinit=ENERGY_THRESHOLD_INIT, valstep=1)

ax_bins = plt.axes([0.2, 0.10, 0.65, 0.03])
slider_bins = Slider(ax_bins, 'Histogram Bins', 20, 200, valinit=BINS_INIT, valstep=1)

def on_slider_change(_=NONE):
    update_plot(slider_short_gate.val, slider_total_gate.val, slider_energy.val, slider_bins.val)

slider_short_gate.on_changed(on_slider_change)
slider_total_gate.on_changed(on_slider_change)
slider_energy.on_changed(on_slider_change)
slider_bins.on_changed(on_slider_change)

update_plot(SHORT_GATE_NS_INIT, TOTAL_GATE_NS_INIT, ENERGY_THRESHOLD_INIT, BINS_INIT)

# Figure 4: Plot Statistical Summary
fig4, ((ax4a, ax4b), (ax4c, ax4d)) = plt.subplots(2, 2, figsize=(12, 10), num="Statistical Analysis - NE213")

# Shows the Amplitude distribution
ax4a.hist(max_amplitudes, bins=50, color='purple', alpha=0.7, edgecolor='black')
ax4a.set_xlabel('Peak Amplitude (a.u.)')
ax4a.set_ylabel('Counts')
ax4a.set_title('Pulse Amplitude Distribution')
ax4a.grid(True, alpha=0.3)

# Shows the Rise time distribution
rise_times_ns = [rt * NS_PER_SAMPLE for rt in rise_times]
ax4b.hist(rise_times_ns, bins=50, color='orange', alpha=0.7, edgecolor='black')
ax4b.set_xlabel('Rise Time (ns)')
ax4b.set_ylabel('Counts')
ax4b.set_title('Pulse Rise Time Distribution\n(10%-90% of peak)')
ax4b.grid(True, alpha=0.3)

# Comparision of Qtot vs Qshort scatter plot, only for valid events
valid_qtot_plot = [qtot[i] for i, v in enumerate(valid_flags) if v]
valid_qshort_plot = [qshort[i] for i, v in enumerate(valid_flags) if v]
valid_psd_plot = [psd[i] for i, v in enumerate(valid_flags) if v]

ax4c.scatter(valid_qtot_plot, valid_qshort_plot, alpha=0.5, s=10, c=valid_psd_plot, cmap='coolwarm')
ax4c.set_xlabel('Qtot (total charge)')
ax4c.set_ylabel('Qshort (short gate charge)')
ax4c.set_title('Short vs Total Charge (Binda Method)\n(Valid events, pile-up rejected)')
ax4c.grid(True, alpha=0.3)

# PSD parameter distribution
ax4d.hist(valid_psd_plot, bins=100, color='green', alpha=0.7, edgecolor='black')
ax4d.set_xlabel('PSD = (Qtot - Qshort) / Qtot')
ax4d.set_ylabel('Counts')
ax4d.set_title('PSD Parameter Distribution (log scale)')
ax4d.set_yscale('log')
ax4d.grid(True, alpha=0.3)

plt.tight_layout()


print("="*70)
print("NE213 NEUTRON DETECTION ANALYSIS - Farnsworth Fusion Reactor")
print("Binda Charge Comparison Method Implementation")
print("="*70)
print(f"Total events captured: {len(waveforms)}")
print(f"Valid events (pile-up rejected): {sum(valid_flags)} ({sum(valid_flags)/len(valid_flags)*100:.1f}%)")
print(f"Pile-up events rejected: {len(valid_flags) - sum(valid_flags)} ({(len(valid_flags) - sum(valid_flags))/len(valid_flags)*100:.1f}%)")
print()
print("Detector Configuration:")
print("  Type: NE213 Liquid Scintillator")
print("  Target reaction: D-D fusion")
print(f"  Sampling rate: {SAMPLE_RATE_GS} GS/s ({NS_PER_SAMPLE} ns/sample)")
print()
print("Pulse Characteristics:")
print(f"  Average amplitude: {np.mean(max_amplitudes):.2f}")
print(f"  Average rise time: {np.mean(rise_times_ns):.2f} ns")
print()
print("Binda Charge Comparison Method (Equation 8):")
print("   PSD = (Qtot - Qshort) / Qtot")
print(f"  Short gate: {SHORT_GATE_NS_INIT} ns")
print(f"  Total gate: {TOTAL_GATE_NS_INIT} ns")
print(f"  PSD range: [{min(valid_psd_plot):.4f}, {max(valid_psd_plot):.4f}]")
print(f"  PSD mean: {np.mean(valid_psd_plot):.4f}")
print(f"  PSD std: {np.std(valid_psd_plot):.4f}")
print()
print("Integrated Charges:")
print(f"  Qtot mean: {np.mean(valid_qtot_plot):.1f}")
print(f"  Qshort mean: {np.mean(valid_qshort_plot):.1f}")
print(f"  Slow component fraction: {(np.mean(valid_qtot_plot) - np.mean(valid_qshort_plot))/np.mean(valid_qtot_plot)*100:.1f}%")
print()
fom, gamma_peak, neutron_peak = calculate_fom(psd, qtot, 10, valid_flags)
print("Discrimination Performance (Binda Section 6.2):")
print(f"  Figure of Merit (FOM): {fom:.4f}")
print("  Reference values:")
print("    - Binda (NE213 at JET): FOM = 2.25-2.35 (VERY GOOD)")
print("    - Baselga (CCM baseline): FOM = 0.94")
print("    - Baselga (best method): FOM = 1.04")
if fom > 2.0:
    print("  [RESULT] VERY GGOOD discrimination (matching JET NE213 performance)")
    neutron_count = sum(1 for psd_val, v in zip(psd, valid_flags) if psd_val > np.median(valid_psd_plot) and v)
    print(f"  [RESULT] Estimated neutron events: {neutron_count} ({neutron_count/sum(valid_flags)*100:.1f}%)")
    print(f"Actual FOM: {fom:.10f}")
elif fom > 1.0:
    print("   [RESULT] GOOD discrimination capability")
    neutron_count = sum(1 for psd_val, v in zip(psd, valid_flags) if psd_val > np.median(valid_psd_plot) and v)
    print(f"  [RESULT] Estimated neutron events: {neutron_count} ({neutron_count/sum(valid_flags)*100:.1f}%)")
    print(f"Actual FOM: {fom:.10f}")
elif fom > 0.7:
    print("  [RESULT] MODERATE discrimination -> Talk to Felix")
    neutron_count = sum(1 for psd_val, v in zip(psd, valid_flags) if psd_val > np.median(valid_psd_plot) and v)
    print(f"  [RESULT] Estimated neutron events: {neutron_count} ({neutron_count/sum(valid_flags)*100:.1f}%)")
    print(f"Actual FOM: {fom:.10f}")
else:
    print("  [RESULT] POOR discrimination -> Talk to Felix")
    print(f"Actual FOM: {fom:.10f}")
if gamma_peak > 0 and neutron_peak > 0:
    print("\nPeak Positions:")
    print(f"  Gamma peak (R): {gamma_peak:.4f}")
    print(f"  Neutron peak (R): {neutron_peak:.4f}")
    print(f"  Separation: {neutron_peak - gamma_peak:.4f}")
    print(f"  Peak ratio: {neutron_peak/gamma_peak:.2f}")
    print()

plt.show()