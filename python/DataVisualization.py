from pickle import NONE
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from scipy import signal
from scipy.optimize import curve_fit
import tkinter as tk
from tkinter import filedialog, ttk
import os

##############################################################
# NOTE: AdrGos: Recommended to run in a virtual environment
# python -m venv ne213_env
# source ne213_env/bin/activate
##############################################################

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

class WaveformLoader:
    """
    Handles loading and classification of waveform files

    Returns:
        _type_: _description_
    """

    @staticmethod
    def load_waveform_file(filename):
        """
        Load waveforms from a text file

        Args:
            filename (str): Path to the waveform data file

        Returns:
            list: List of waveforms, each waveform is a list of float values
        """
        with open(filename, 'r') as f:
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

        print(f"Loaded {len(waveforms)} waveforms from file: {filename}")
        return waveforms

    @staticmethod
    def classify_file_type(filename):
        """
        Classify file as 'mit' (with fusor) or 'ohne' (background)

        Args:
            filename (str): Path to the waveform data file

        Returns:
            str: 'mit' if file contains fusor data, 'ohne' if background data
        """
        basename = os.path.basename(filename).lower()
        if 'mit' in basename:
            return 'mit'
        elif 'ohne' in basename:
            return 'ohne'
        else:
            return None

class FileSelectionGUI:
    """
    GUI for file selection and parameter configuration
    """

    def __init__(self):
        """
        Initialize the GUI components
        """
        self.root = tk.Tk()
        self.root.title("NE213 Analyzer - Configuration")
        self.root.geometry("600x450")

        self.selected_file = None
        self.short_gate_ns = 25
        self.total_gate_ns = 90
        self.file_type = None
        self.preview_processor = None

        self._create_widgets()

    def _create_widgets(self):
        """
        Create and layout GUI widgets
        """
        main_frame = ttk.Frame(self.root, padding="20")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        title_label = ttk.Label(main_frame, text="NE213 Waveform Analyzer", font=('Arial', 16, 'bold'))
        title_label.grid(row=0, column=0, columnspan=3, pady=(0, 20))

        file_frame = ttk.LabelFrame(main_frame, text="Step 1: File Selection", padding="10")
        file_frame.grid(row=1, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(0, 20))

        self.file_label = ttk.Label(file_frame, text="No file selected", foreground="gray")
        self.file_label.grid(row=0, column=0, sticky=tk.W, padx=(0, 10))

        browse_btn = ttk.Button(file_frame, text="Browse...", command=self._browse_file)
        browse_btn.grid(row=0, column=1)

        self.type_label = ttk.Label(file_frame, text="", foreground="blue")
        self.type_label.grid(row=1, column=0, columnspan=2, sticky=tk.W, pady=(5, 0))

        preview_frame = ttk.Frame(main_frame)
        preview_frame.grid(row=2, column=0, columnspan=3, pady=(0, 10))

        self.preview_btn = ttk.Button(preview_frame, text="Step 2: Preview Waveforms", command=self._preview_waveforms, state=tk.DISABLED)
        self.preview_btn.pack()

        ttk.Label(preview_frame, text="(Load and view waveforms to determine gate parameters)", foreground="gray", font=('Arial', 9)).pack(pady=(5, 0))

        params_frame = ttk.LabelFrame(main_frame, text="Step 3: Integration Gate Parameters", padding="10")
        params_frame.grid(row=3, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(0, 20))

        ttk.Label(params_frame, text="Short Gate (ns):").grid(row=0, column=0, sticky=tk.W, pady=5)
        self.short_gate_var = tk.StringVar(value="25")
        short_gate_entry = ttk.Entry(params_frame, textvariable=self.short_gate_var, width=10)
        short_gate_entry.grid(row=0, column=1, padx=(10, 20), pady=5)
        ttk.Label(params_frame, text="(10-100 ns)", foreground="gray").grid(row=0, column=2, sticky=tk.W)

        ttk.Label(params_frame, text="Total Gate (ns):").grid(row=1, column=0, sticky=tk.W, pady=5)
        self.total_gate_var = tk.StringVar(value="90")
        total_gate_entry = ttk.Entry(params_frame, textvariable=self.total_gate_var, width=10)
        total_gate_entry.grid(row=1, column=1, padx=(10, 20), pady=5)
        ttk.Label(params_frame, text="(50-300 ns)", foreground="gray").grid(row=1, column=2, sticky=tk.W)

        info_frame = ttk.Frame(main_frame)
        info_frame.grid(row=4, column=0, columnspan=3, pady=(0, 20))

        info_text = "Binda recommended: Short=25ns, Total=90ns"
        ttk.Label(info_frame, text=info_text, foreground="gray", justify=tk.CENTER).pack()

        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=5, column=0, columnspan=3)

        self.start_btn = ttk.Button(button_frame, text="Step 4: Start Full Analysis", command=self._start_analysis, state=tk.DISABLED)
        self.start_btn.pack(side=tk.LEFT, padx=5)

        cancel_btn = ttk.Button(button_frame, text="Cancel", command=self.root.quit)
        cancel_btn.pack(side=tk.LEFT, padx=5)

    def _browse_file(self):
        """
        Open file dialog to select waveform file
        """
        file_path = filedialog.askopenfilename(
            title="Select NE213 waveform file",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )

        if file_path:
            self.selected_file = file_path
            filename = os.path.basename(file_path)
            self.file_label.config(text=filename, foreground="black")

            self.file_type = WaveformLoader.classify_file_type(file_path)
            if self.file_type == 'mit':
                self.type_label.config(text="File type: WITH fusor (mit)", foreground="green")
            elif self.file_type == 'ohne':
                self.type_label.config(text="File type: BACKGROUND (ohne)", foreground="orange")
            else:
                self.type_label.config(text="File type: UNKNOWN - classification failed", foreground="red")

            self.preview_btn.config(state=tk.NORMAL)

    def _preview_waveforms(self):
        if not self.selected_file:
            return

        self.preview_processor = WaveformProcessor()
        waveforms = WaveformLoader.load_waveform_file(self.selected_file)
        self.preview_processor.waveforms = waveforms
        self.preview_processor._clean_waveforms()
        self.preview_processor._baseline_correction()

        fig_preview = plt.figure("Waveform Preview - Determine Gate Parameters", figsize=(14, 6))
        ax1 = fig_preview.add_subplot(121)
        time_axis = np.arange(350) * NS_PER_SAMPLE

        for i, waveform in enumerate(self.preview_processor.waveforms_baseline_corrected[:50]):
            ax1.plot(time_axis[:len(waveform)], waveform, alpha=0.3, linewidth=0.5)

        ax1.set_xlabel('Time (ns)')
        ax1.set_ylabel('Amplitude (baseline corrected)')
        ax1.set_title('Raw Waveforms Preview (first 50 events)\nUse this to determine Short and Total gate positions')
        ax1.grid(True, alpha=0.3)
        ax1.axvline(x=25, color='green', linestyle='--', alpha=0.5, label='Suggested Short gate (25 ns)')
        ax1.axvline(x=90, color='red', linestyle='--', alpha=0.5, label='Suggested Total gate (90 ns)')
        ax1.legend()

        ax2 = fig_preview.add_subplot(122)
        if len(self.preview_processor.waveforms_baseline_corrected) > 0:
            avg_waveform = np.mean(self.preview_processor.waveforms_baseline_corrected, axis=0)
            ax2.plot(time_axis[:len(avg_waveform)], avg_waveform, 'b-', linewidth=2, label='Average waveform')
            ax2.set_xlabel('Time (ns)')
            ax2.set_ylabel('Amplitude')
            ax2.set_title('Average Waveform\nShort gate captures fast component, Total gate captures slow tail')
            ax2.legend()
            ax2.grid(True, alpha=0.3)
            ax2.axvline(x=25, color='green', linestyle='--', alpha=0.5)
            ax2.axvline(x=90, color='red', linestyle='--', alpha=0.5)
            ax2.set_yscale('log')
            ax2.set_ylim(bottom=0.1)

        plt.tight_layout()
        plt.show(block=False)

        self.start_btn.config(state=tk.NORMAL)

    def _start_analysis(self):
        """
        Validate parameters and close GUI to start analysis
        """
        try:
            self.short_gate_ns = float(self.short_gate_var.get())
            self.total_gate_ns = float(self.total_gate_var.get())

            if not (10 <= self.short_gate_ns <= 100):
                tk.messagebox.showerror("Invalid Input", "Short gate must be between 10 and 100 ns")
                return

            if not (50 <= self.total_gate_ns <= 300):
                tk.messagebox.showerror("Invalid Input", "Total gate must be between 50 and 300 ns")
                return

            if self.short_gate_ns >= self.total_gate_ns:
                tk.messagebox.showerror("Invalid Input", "Short gate must be less than total gate")
                return

            plt.close('all')
            self.root.quit()

        except ValueError:
            tk.messagebox.showerror("Invalid Input", "Please enter valid numbers for gate parameters")

    def show(self):
        """
        Show the GUI and return selected parameters
        :return: A tuple (selected_file, short_gate_ns, total_gate_ns)
        """
        self.root.mainloop()
        return self.selected_file, self.short_gate_ns, self.total_gate_ns

class WaveformProcessor:
    """
    Processes and analyzes waveform data
    """

    def __init__(self, short_gate_ns=25, total_gate_ns=90):
        """
        Initalizes all necessary variables
        """
        self.waveforms = []
        self.waveforms_clean = []
        self.waveforms_baseline_corrected = []
        self.baselines = []
        self.max_amplitudes = []
        self.rise_times = []
        self.max_indices = []
        self.psd_values = []
        self.qtot_values = []
        self.qshort_values = []
        self.valid_flags = []
        self.file_type = None
        self.short_gate_ns = short_gate_ns
        self.total_gate_ns = total_gate_ns

    def load_and_process(self, filename):
        """
        Load and process waveforms from file

        Args:
            filename (str): Path to the waveform data file
        """

        # Load waveforms and classify file type...
        self.file_type = WaveformLoader.classify_file_type(filename)
        self.waveforms = WaveformLoader.load_waveform_file(filename)

        # Do signal preprocessing...
        self._clean_waveforms()
        self._baseline_correction()

        if self.file_type == 'mit':
            # Only calculate PSD for 'mit' files (with fusor running)
            self._calculate_psd_parameters()
        else:
            # For 'ohne' files (background), we don't need PSD calculations -> because no neutrons expected in background measurements
            print("Background file detected - skipping PSD calculations")

    def _clean_waveforms(self):
        """
        Clean and invert waveforms
        """
        self.waveforms_clean = []
        for waveform_orig in self.waveforms:
            waveform_inverted = [-x for x in waveform_orig]
            self.waveforms_clean.append(waveform_inverted)

    def _baseline_correction(self):
        """
        Apply baseline correction to waveforms with adaptive noise filtering
        """
        self.waveforms_baseline_corrected = []
        self.baselines = []
        self.max_amplitudes = []
        self.rise_times = []
        self.max_indices = []
        self.rejected_waveforms = []

        baseline_stds = []
        for waveform in self.waveforms_clean:
            baseline_region = waveform[:50]
            baseline_stds.append(np.std(baseline_region))

        # Adaptive threshold: reject waveforms with baseline noise > 2x median std
        median_baseline_std = np.median(baseline_stds)
        noise_threshold = median_baseline_std * 3.0  # 3-sigma threshold

        print(f"Baseline noise analysis:")
        print(f"  Median baseline std: {median_baseline_std:.2f}")
        print(f"  Adaptive noise threshold: {noise_threshold:.2f}")

        for i, waveform in enumerate(self.waveforms_clean):
            # Baseline correction using first 50 samples
            baseline_region = waveform[:50]
            baseline = np.mean(baseline_region)
            baseline_std = np.std(baseline_region)

            self.baselines.append(baseline)
            waveform_baseline_corrected = [(x - baseline) for x in waveform]

            # Adaptive noise rejection based on baseline statistics
            if baseline_std > noise_threshold:
                self.rejected_waveforms.append({
                    'index': i,
                    'reason': 'high_baseline_noise',
                    'baseline_std': baseline_std
                })
                continue

            # Check for unphysical negative spikes in pre-trigger region
            # These indicate electrical noise, not real scintillation
            first_hunderd_samples = waveform_baseline_corrected[:100]
            min_pretrigger = min(first_hunderd_samples)

            # Adaptive threshold based on baseline noise level
            if min_pretrigger < -(6 * baseline_std + 100):  # 5-sigma + safety margin
                self.rejected_waveforms.append({
                    'index': i,
                    'reason': 'negative_spike',
                    'min_value': min_pretrigger,
                    'threshold': -(6 * baseline_std + 100)
                })
                continue

            self.waveforms_baseline_corrected.append(waveform_baseline_corrected)

            # Find peak amplitude and index
            max_amp = max(waveform_baseline_corrected)
            self.max_amplitudes.append(max_amp)
            max_idx = np.argmax(waveform_baseline_corrected)
            self.max_indices.append(max_idx)

            # Calculate rise time, we'll try 10% to 90% of peak
            threshold_10 = 0.1 * max_amp
            threshold_90 = 0.9 * max_amp
            idx_10 = next((i for i, v in enumerate(waveform_baseline_corrected) if v > threshold_10), None)
            idx_90 = next((i for i, v in enumerate(waveform_baseline_corrected) if v > threshold_90), None)
            if idx_10 is not None and idx_90 is not None:
                self.rise_times.append(idx_90 - idx_10)
            else:
                self.rise_times.append(0)

        print(f"  Waveforms rejected: {len(self.rejected_waveforms)} / {len(self.waveforms_clean)} ({len(self.rejected_waveforms)/len(self.waveforms_clean)*100:.1f}%)")
        if self.rejected_waveforms:
            noise_rejects = sum(1 for r in self.rejected_waveforms if r['reason'] == 'high_baseline_noise')
            spike_rejects = sum(1 for r in self.rejected_waveforms if r['reason'] == 'negative_spike')
            print(f"    - High baseline noise: {noise_rejects}")
            print(f"    - Negative spikes: {spike_rejects}")

    def detect_pileup(self, waveform, max_idx, qtot, psd_factor, threshold_factor=0.25):
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
        if psd_factor > 0.4:  #TODO AdrGos: Verify with Felix if this threshold is valid for NE213 or better to use 0.4-0.45 to reject more pile-up
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

    def _calculate_psd_parameters(self, short_gate_ns=None, total_gate_ns=None, pre_trigger_ns=10):
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
            short_gate_ns: Duration of short integration gate (ns)
            total_gate_ns: Duration of total integration gate (ns)
            pre_trigger_ns: Time before peak to start integration (ns)
        """
        if short_gate_ns is None:
            short_gate_ns = self.short_gate_ns
        if total_gate_ns is None:
            total_gate_ns = self.total_gate_ns

        short_gate_samples = ns_to_samples(short_gate_ns)
        total_gate_samples = ns_to_samples(total_gate_ns)
        pre_trigger_samples = ns_to_samples(pre_trigger_ns)

        self.psd_values = []
        self.qtot_values = []
        self.qshort_values = []
        self.valid_flags = []

        for i, (waveform, max_idx) in enumerate(zip(self.waveforms_baseline_corrected, self.max_indices)):
            if max_idx is None or max_idx >= len(waveform) - 10:
                self.psd_values.append(0)
                self.qtot_values.append(0)
                self.qshort_values.append(0)
                self.valid_flags.append(False)
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
            is_pileup = self.detect_pileup(waveform, max_idx, qtot, psd)

            self.psd_values.append(psd)
            self.qtot_values.append(qtot)
            self.qshort_values.append(qshort)
            self.valid_flags.append(not is_pileup)

        return self.psd_values, self.qtot_values, self.qshort_values, self.valid_flags

    def calculate_fom(self, psd_values, qtot_values, valid_flags, energy_threshold=10):
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
            energy_threshold: Minimum energy for discrimination

        Returns:
            fom: Figure of merit
            gamma_peak: Peak position for gamma events
            neutron_peak: Peak position for neutron events
        """
        if self.file_type != 'mit':
            print("FOM calculation only applicable for 'mit' files (with fusor)")
            return 0, 0, 0

        # Filter by energy threshold and validity
        valid_idx = [i for i, (e, v) in enumerate(zip(qtot_values, valid_flags))
                     if e > energy_threshold and v]

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

class NE213Analyzer:
    """
    Main analyzer class for NE213 data visualization and analysis
    """

    def __init__(self, short_gate_ns=25, total_gate_ns=90):
        """
        Initialize the NE213Analyzer with a waveform processor and figure list
        """
        self.processor = WaveformProcessor(short_gate_ns, total_gate_ns)
        self.figures = []
        self.short_gate_ns = short_gate_ns
        self.total_gate_ns = total_gate_ns

    def select_and_load_file(self):
        """
        Open enhanced GUI for file selection and parameter configuration

        Returns:
            bool: True if file loaded successfully, False otherwise
        """
        gui = FileSelectionGUI()
        file_path, self.short_gate_ns, self.total_gate_ns = gui.show()
        gui.root.destroy()

        if file_path:
            self.processor = WaveformProcessor(self.short_gate_ns, self.total_gate_ns)
            self.processor.load_and_process(file_path)
            return True
        return False

    def plot_raw_waveforms(self):
        """
        Figure 1: Plot Raw Waveforms
        """
        fig1 = plt.figure("Raw Waveforms - NE213 Analysis", figsize=(14, 5))
        ax1 = fig1.add_subplot(121)
        time_axis = np.arange(350) * NS_PER_SAMPLE  # Time axis in nanoseconds

        # Plot all waveforms...
        for i, waveform in enumerate(self.processor.waveforms_baseline_corrected):
            ax1.plot(time_axis[:len(waveform)], waveform, alpha=0.3, linewidth=0.5)

        ax1.set_xlabel('Time (ns)')
        ax1.set_ylabel('Amplitude (baseline corrected)')
        ax1.set_title('NE213 Scintillation Pulses (D-D Fusion)')
        ax1.grid(True, alpha=0.3)
        ax1.axvline(x=self.short_gate_ns, color='green', linestyle='--', alpha=0.5, label=f'Short gate end ({self.short_gate_ns} ns)')
        ax1.axvline(x=self.total_gate_ns, color='red', linestyle='--', alpha=0.5, label=f'Total gate end ({self.total_gate_ns} ns)')
        ax1.legend()

        # Plot Average pulse shapes comparison from binda...
        ax2 = fig1.add_subplot(122)
        if self.processor.file_type == 'mit':
            sorted_by_psd = sorted(zip(self.processor.psd_values, self.processor.waveforms_baseline_corrected, self.processor.valid_flags), key=lambda x: x[0])
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
                    ax2.axvline(x=self.short_gate_ns, color='green', linestyle='--', alpha=0.3)
                    ax2.axvline(x=self.total_gate_ns, color='red', linestyle='--', alpha=0.3)
        else:
            # For background files, just plot average waveform...
            if self.processor.waveforms_baseline_corrected:
                avg_waveform = np.mean(self.processor.waveforms_baseline_corrected, axis=0)
                ax2.plot(time_axis[:len(avg_waveform)], avg_waveform, 'g-', label='Background Average', linewidth=2)
                ax2.set_xlabel('Time (ns)')
                ax2.set_ylabel('Amplitude')
                ax2.set_title('Average Background Pulse Shape')
                ax2.legend()
                ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        self.figures.append(fig1)

    def plot_psd_analysis(self):
        """
        Figure 2: creates 2D PSD Plot (only for mit files
        """
        if self.processor.file_type != 'mit':
            return

        fig2, ax_psd = plt.subplots(figsize=(10, 8), num="PSD Analysis: NE213 Neutron-Gamma Discrimination")
        # Filter for valid events
        valid_qtot = [self.processor.qtot_values[i] for i, v in enumerate(self.processor.valid_flags) if v]
        valid_psd = [self.processor.psd_values[i] for i, v in enumerate(self.processor.valid_flags) if v]

        scatter = ax_psd.scatter(valid_qtot, valid_psd, c=valid_psd, cmap='coolwarm', alpha=0.6, s=20, edgecolors='black', linewidth=0.5)
        ax_psd.set_xlabel('Qtot - Total Charge (Energy Proxy)', fontsize=12)
        ax_psd.set_ylabel('PSD = (Qtot - Qshort) / Qtot', fontsize=12)
        ax_psd.set_title(f'NE213 Pulse Shape Discrimination (Binda Method)\n' +
                         f'Short={self.short_gate_ns}ns, Total={self.total_gate_ns}ns\n' +
                         'D-D Fusion\n' +
                         '(Neutrons: high PSD, Gammas: low PSD)', fontsize=14)
        ax_psd.grid(True, alpha=0.3)
        cbar = plt.colorbar(scatter, ax=ax_psd)
        cbar.set_label('PSD Parameter', fontsize=10)

        # Add discrimination line
        if len(valid_psd) > 0:
            median_psd = np.median(valid_psd)
            ax_psd.axhline(y=median_psd, color='green', linestyle='--', linewidth=2, label=f'Discrimination threshold (median PSD={median_psd:.3f})')
            ax_psd.legend()

        self.figures.append(fig2)

    def plot_interactive_histogram(self):
        """
        Figure 3: Interactive Histogram with Parameter Control (only for mit files)
        """
        if self.processor.file_type != 'mit':
            return

        fig3 = plt.figure("Interactive PSD Histogram - NE213 Optimization", figsize=(12, 8))
        plt.subplots_adjust(bottom=0.35)
        ax_hist = fig3.add_subplot(111)

        # Initial parameters in NANOSECONDS
        SHORT_GATE_NS_INIT = self.short_gate_ns
        TOTAL_GATE_NS_INIT = self.total_gate_ns
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
            # Create a temporary processor instance to avoid overwriting original data
            temp_processor = WaveformProcessor()
            temp_processor.waveforms_baseline_corrected = self.processor.waveforms_baseline_corrected
            temp_processor.max_indices = self.processor.max_indices
            temp_processor.file_type = self.processor.file_type

            psd_new, qtot_new, _, valid_new = temp_processor._calculate_psd_parameters(
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
                ax_hist.hist(filtered_psd, bins=int(nbins), color='blue', alpha=0.7, edgecolor='black')
                fom, gamma_pk, neutron_pk = temp_processor.calculate_fom(psd_new, qtot_new, valid_new, energy_threshold)

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
        ax_short_gate = fig3.add_axes([0.2, 0.25, 0.65, 0.03])
        slider_short_gate = Slider(ax_short_gate, 'Short Gate (ns)', 10, 100, valinit=SHORT_GATE_NS_INIT, valstep=0.5)

        ax_total_gate = fig3.add_axes([0.2, 0.20, 0.65, 0.03])
        slider_total_gate = Slider(ax_total_gate, 'Total Gate (ns)', 50, 300, valinit=TOTAL_GATE_NS_INIT, valstep=1)

        ax_energy = fig3.add_axes([0.2, 0.15, 0.65, 0.03])
        slider_energy = Slider(ax_energy, 'Qtot Threshold (energy)', 0, 100, valinit=ENERGY_THRESHOLD_INIT, valstep=1)

        ax_bins = fig3.add_axes([0.2, 0.10, 0.65, 0.03])
        slider_bins = Slider(ax_bins, 'Histogram Bins', 20, 200, valinit=BINS_INIT, valstep=1)


        def on_slider_change(_=None):
            """
            Update the plot when any slider is changed.
            """
            update_plot(
                slider_short_gate.val,
                slider_total_gate.val,
                slider_energy.val,
                slider_bins.val
            )

        slider_short_gate.on_changed(on_slider_change)
        slider_total_gate.on_changed(on_slider_change)
        slider_energy.on_changed(on_slider_change)
        slider_bins.on_changed(on_slider_change)

        # Initial plot
        update_plot(SHORT_GATE_NS_INIT, TOTAL_GATE_NS_INIT, ENERGY_THRESHOLD_INIT, BINS_INIT)
        fig3.sliders = [slider_short_gate, slider_total_gate, slider_energy, slider_bins]
        self.figures.append(fig3)

    def plot_statistical_summary(self):
        """
        Figure 4: Plot Statistical Summary
        """
        fig4, ((ax4a, ax4b), (ax4c, ax4d)) = plt.subplots(2, 2, figsize=(12, 10), num="Statistical Analysis - NE213")

        # Shows the Amplitude distribution
        ax4a.hist(self.processor.max_amplitudes, bins=50, color='purple', alpha=0.7, edgecolor='black')
        ax4a.set_xlabel('Peak Amplitude (a.u.)')
        ax4a.set_ylabel('Counts')
        ax4a.set_title('Pulse Amplitude Distribution')
        ax4a.grid(True, alpha=0.3)

        # Shows the Rise time distribution
        rise_times_ns = [rt * NS_PER_SAMPLE for rt in self.processor.rise_times]
        ax4b.hist(rise_times_ns, bins=50, color='orange', alpha=0.7, edgecolor='black')
        ax4b.set_xlabel('Rise Time (ns)')
        ax4b.set_ylabel('Counts')
        ax4b.set_title('Pulse Rise Time Distribution\n(10%-90% of peak)')
        ax4b.grid(True, alpha=0.3)

        if self.processor.file_type == 'mit':
            # Comparision of Qtot vs Qshort scatter plot, only for valid events
            valid_qtot_plot = [self.processor.qtot_values[i] for i, v in enumerate(self.processor.valid_flags) if v]
            valid_qshort_plot = [self.processor.qshort_values[i] for i, v in enumerate(self.processor.valid_flags) if v]
            valid_psd_plot = [self.processor.psd_values[i] for i, v in enumerate(self.processor.valid_flags) if v]

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
        else:
            # For background files
            ax4c.text(0.5, 0.5, 'Background Measurement\nNo PSD Analysis Required',
                      horizontalalignment='center', verticalalignment='center',
                      transform=ax4c.transAxes, fontsize=12)
            ax4c.set_title('Background Measurement')
            ax4c.axis('off')

            ax4d.text(0.5, 0.5, 'Use this data for\nnoise/background characterization',
                      horizontalalignment='center', verticalalignment='center',
                      transform=ax4d.transAxes, fontsize=12)
            ax4d.set_title('Background Reference')
            ax4d.axis('off')

        plt.tight_layout()
        self.figures.append(fig4)

    def print_summary(self):
        """
        Print analysis summary
        """
        print("="*70)
        print("NE213 NEUTRON DETECTION ANALYSIS - Farnsworth Fusion Reactor")
        print(f"File type: {'WITH fusor (mit)' if self.processor.file_type == 'mit' else 'BACKGROUND (ohne)'}")
        print("Binda Charge Comparison Method Implementation")
        print("="*70)
        print(f"Total events captured: {len(self.processor.waveforms)}")

        if self.processor.file_type == 'mit':
            print(f"Valid events (pile-up rejected): {sum(self.processor.valid_flags)} ({sum(self.processor.valid_flags)/len(self.processor.valid_flags)*100:.1f}%)")
            print(f"Pile-up events rejected: {len(self.processor.valid_flags) - sum(self.processor.valid_flags)} ({(len(self.processor.valid_flags) - sum(self.processor.valid_flags))/len(self.processor.valid_flags)*100:.1f}%)")
        else:
            print("Background measurement - no pile-up rejection applied")

        print()
        print("Detector Configuration:")
        print("  Type: NE213 Liquid Scintillator")
        print("  Target reaction: D-D fusion")
        print(f"  Sampling rate: {SAMPLE_RATE_GS} GS/s ({NS_PER_SAMPLE} ns/sample)")
        print()
        print("Pulse Characteristics:")
        print(f"  Average amplitude: {np.mean(self.processor.max_amplitudes):.2f}")
        print(f"  Average rise time: {np.mean([rt * NS_PER_SAMPLE for rt in self.processor.rise_times]):.2f} ns")
        print()

        if self.processor.file_type == 'mit':
            print("Binda Charge Comparison Method (Equation 8):")
            print("   PSD = (Qtot - Qshort) / Qtot")
            print(f"  Short gate: {self.short_gate_ns} ns")
            print(f"  Total gate: {self.total_gate_ns} ns")

            valid_psd_plot = [self.processor.psd_values[i] for i, v in enumerate(self.processor.valid_flags) if v]
            if valid_psd_plot:
                print(f"  PSD range: [{min(valid_psd_plot):.4f}, {max(valid_psd_plot):.4f}]")
                print(f"  PSD mean: {np.mean(valid_psd_plot):.4f}")
                print(f"  PSD std: {np.std(valid_psd_plot):.4f}")
            print()

            valid_qtot_plot = [self.processor.qtot_values[i] for i, v in enumerate(self.processor.valid_flags) if v]
            valid_qshort_plot = [self.processor.qshort_values[i] for i, v in enumerate(self.processor.valid_flags) if v]
            print("Integrated Charges:")
            if valid_qtot_plot and valid_qshort_plot:
                print(f"  Qtot mean: {np.mean(valid_qtot_plot):.1f}")
                print(f"  Qshort mean: {np.mean(valid_qshort_plot):.1f}")
                print(f"  Slow component fraction: {(np.mean(valid_qtot_plot) - np.mean(valid_qshort_plot))/np.mean(valid_qtot_plot)*100:.1f}%")
            print()

            fom, gamma_peak, neutron_peak = self.processor.calculate_fom(
                self.processor.psd_values,
                self.processor.qtot_values,
                self.processor.valid_flags,
                10
            )
            print("Discrimination Performance (Binda Section 6.2):")
            print(f"  Figure of Merit (FOM): {fom:.4f}")
            print("  Reference values:")
            print("    - Binda (NE213 at JET): FOM = 2.25-2.35 (VERY GOOD)")
            print("    - Baselga (CCM baseline): FOM = 0.94")
            print("    - Baselga (best method): FOM = 1.04")

            if fom > 2.0:
                print("  [RESULT] VERY GGOOD discrimination (matching JET NE213 performance)")
                if valid_psd_plot:
                    median_psd = np.median(valid_psd_plot)
                    neutron_count = sum(1 for psd_val, v in zip(self.processor.psd_values, self.processor.valid_flags) if psd_val > median_psd and v)
                    print(f"  [RESULT] Estimated neutron events: {neutron_count} ({neutron_count/sum(self.processor.valid_flags)*100:.1f}%)")
                print(f"Actual FOM: {fom:.10f}")
            elif fom > 1.0:
                print("   [RESULT] GOOD discrimination capability")
                if valid_psd_plot:
                    median_psd = np.median(valid_psd_plot)
                    neutron_count = sum(1 for psd_val, v in zip(self.processor.psd_values, self.processor.valid_flags) if psd_val > median_psd and v)
                    print(f"  [RESULT] Estimated neutron events: {neutron_count} ({neutron_count/sum(self.processor.valid_flags)*100:.1f}%)")
                print(f"Actual FOM: {fom:.10f}")
            elif fom > 0.7:
                print("  [RESULT] MODERATE discrimination -> Talk to Felix")
                if valid_psd_plot:
                    median_psd = np.median(valid_psd_plot)
                    neutron_count = sum(1 for psd_val, v in zip(self.processor.psd_values, self.processor.valid_flags) if psd_val > median_psd and v)
                    print(f"  [RESULT] Estimated neutron events: {neutron_count} ({neutron_count/sum(self.processor.valid_flags)*100:.1f}%)")
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
        else:
            print("Background measurement summary:")
            print("  This data represents background noise (cosmic rays, natural radioactivity, etc.)")
            print("  Use this for noise characterization and background subtraction")
            print("  No neutron/gamma discrimination analysis performed on background data")

    def run_analysis(self):
        """
        Run complete analysis pipeline
        """
        if not self.select_and_load_file():
            print("No file selected. Exiting.")
            return

        self.plot_raw_waveforms()
        self.plot_psd_analysis()
        self.plot_interactive_histogram()
        self.plot_statistical_summary()
        self.print_summary()

        plt.show()

def main():
    """
    Main function to run the NE213 analyzer
    """
    analyzer = NE213Analyzer()
    analyzer.run_analysis()

if __name__ == "__main__":
    main()