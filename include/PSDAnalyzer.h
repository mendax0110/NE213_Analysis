#pragma once

#include "WaveformData.h"
#include <vector>

/// @brief Namespace for NE213 detector data analysis \namespace ne213
namespace ne213
{
    /**
     * @brief Pulse Shape Discrimination analyzer for NE213 detector
     *
     * Calculates PSD parameters using Binda's Charge Comparison Method (Equation 8):
     *   PSD = (Qtot - Qshort) / Qtot
     *
     * Where:
     * - Qtot: Total charge integrated over full gate (captures entire pulse)
     * - Qshort: Charge integrated over short gate (captures mainly fast component)
     * - The slow component fraction indicates neutron vs gamma
     *
     * Also implements pile-up detection to reject corrupted events.
     */
    class PSDAnalyzer
    {
    public:
        /**
         * @brief Construct a PSD analyzer with specified gate settings
         *
         * Default values are optimized for NE213 with D-D fusion neutrons (2.45 MeV):
         * - Short gate: 25 ns (captures fast scintillation component)
         * - Total gate: 150 ns (captures full pulse including slow component)
         * - Pre-trigger: 10 ns (start integration before peak)
         *
         * @param short_gate_ns Short integration gate duration in nanoseconds
         * @param total_gate_ns Total integration gate duration in nanoseconds
         * @param pre_trigger_ns Pre-trigger time before peak to start integration
         */
        explicit PSDAnalyzer(double short_gate_ns = 25.0, double total_gate_ns = 150.0, double pre_trigger_ns = 10.0);

        /**
         * @brief Calculate PSD parameters for a single waveform
         *
         * Computes Qtot, Qshort, and PSD value for the given waveform.
         * Also performs pile-up detection to flag invalid events.
         *
         * @param waveform Processed waveform data (baseline-corrected)
         * @return Calculated PSD parameters including validity flag
         */
        [[nodiscard]] PSDParameters calculate_parameters(const WaveformData& waveform) const;

        /**
         * @brief Calculate Figure of Merit for neutron-gamma discrimination
         *
         * Implements Binda Equation (9):
         *   FOM = separation / (FWHM_gamma + FWHM_neutron)
         *
         * Reference values (from literature):
         * - Binda (NE213 at JET): FOM = 2.25-2.35 (VERY GOOD)
         * - Baselga (CCM baseline): FOM = 0.94
         * - Baselga (best method): FOM = 1.04
         *
         * @param psd_params Vector of PSD parameters from all waveforms
         * @param energy_threshold Minimum Qtot for event inclusion
         * @return FOM result with peak positions
         */
        static FOMResult calculate_fom(const std::vector<PSDParameters>& psd_params, double energy_threshold = 100);

        /**
         * @brief Create histogram of PSD values above energy threshold
         *
         * Generates a histogram of valid PSD values for events above the
         * energy threshold. Used for FOM calculation and visualization.
         *
         * @param psd_params Vector of PSD parameters
         * @param energy_threshold Minimum Qtot for histogram inclusion
         * @param n_bins Number of histogram bins
         * @return Generated histogram structure
         */
        static Histogram create_histogram(const std::vector<PSDParameters>& psd_params, double energy_threshold, int n_bins = 100);

        /**
         * @brief Set short integration gate duration
         * @param ns Gate duration in nanoseconds (typical: 20-40 ns)
         */
        void set_short_gate(const double ns) { short_gate_ns_ = ns; }

        /**
         * @brief Set total integration gate duration
         * @param ns Gate duration in nanoseconds (typical: 100-200 ns)
         */
        void set_total_gate(const double ns) { total_gate_ns_ = ns; }

        /**
         * @brief Set pre-trigger integration start time
         * @param ns Pre-trigger time in nanoseconds (typical: 5-15 ns)
         */
        void set_pre_trigger(const double ns) { pre_trigger_ns_ = ns; }

        /**
         * @brief Get current short gate duration
         * @return Short gate duration in nanoseconds
         */
        [[nodiscard]] double get_short_gate() const { return short_gate_ns_; }

        /**
         * @brief Get current total gate duration
         * @return Total gate duration in nanoseconds
         */
        [[nodiscard]] double get_total_gate() const { return total_gate_ns_; }

        /**
         * @brief Get current pre-trigger duration
         * @return Pre-trigger duration in nanoseconds
         */
        [[nodiscard]] double get_pre_trigger() const { return pre_trigger_ns_; }

    private:
        double short_gate_ns_;
        double total_gate_ns_;
        double pre_trigger_ns_;

        /**
         * @brief Detect pile-up events using Binda's criteria
         *
         * Implements pile-up detection from Binda Section 8.1.2:
         * 1. Abnormally high PSD (> 0.5) indicates pile-up
         * 2. Multiple peaks in rising edge indicate overlapping pulses
         * 3. Abnormal decay (sudden increases) indicates pile-up
         *
         * @param waveform Baseline-corrected waveform data
         * @param max_idx Index of peak amplitude
         * @param qtot Total integrated charge
         * @param psd_value Calculated PSD parameter
         * @return true if pile-up detected (event should be rejected)
         */
        [[nodiscard]] bool detect_pileup(const std::vector<double>& waveform, int max_idx, double qtot, double psd_value) const;

        /**
         * @brief Find peaks in histogram for FOM calculation
         *
         * Identifies gamma and neutron peaks in PSD histogram using
         * height threshold and prominence criteria.
         *
         * @param counts Histogram bin counts
         * @param min_height Minimum relative height threshold (fraction of max)
         * @return Indices of detected peaks sorted by position
         */
        static std::vector<int> find_peaks(const std::vector<int>& counts, double min_height = 0.15);

        /**
         * @brief Calculate Full Width at Half Maximum of a histogram peak
         *
         * Measures FWHM for FOM calculation (Binda Equation 9).
         *
         * @param counts Histogram bin counts
         * @param bin_centers Center values of histogram bins
         * @param peak_idx Index of the peak in the histogram
         * @return FWHM value in PSD units
         */
        static double calculate_fwhm(const std::vector<int>& counts, const std::vector<double>& bin_centers, int peak_idx);
    };
}