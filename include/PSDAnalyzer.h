#pragma once

#include "WaveformData.h"
#include <vector>

/// @brief Namespace for NE213 detector data analysis \namespace ne213
namespace ne213
{
    /// @brief Class for Pulse Shape Discrimination (PSD) analysis \class PSDAnalyzer
    class PSDAnalyzer
    {
    public:
        /**
         * @brief Constructor to initialize PSDAnalyzer with gate settings
         * @param short_gate_ns Short integration gate duration in nanoseconds
         * @param total_gate_ns Total integration gate duration in nanoseconds
         * @param pre_trigger_ns Pre-trigger duration in nanoseconds
         */
        PSDAnalyzer(double short_gate_ns = 25.0, double total_gate_ns = 150.0, double pre_trigger_ns = 10.0);

        /**
         * @brief Calculate PSD parameters for a given waveform
         * @param waveform Processed waveform data
         * @return Calculated PSD parameters
         */
        PSDParameters calculate_parameters(const WaveformData& waveform) const;

        /**
         * @brief Calculate Figure of Merit (FOM) from PSD parameters
         * @param psd_params Vector of PSD parameters
         * @param energy_threshold Energy threshold for FOM calculation
         * @return Calculated FOM result
         */
        static FOMResult calculate_fom(const std::vector<PSDParameters>& psd_params, double energy_threshold = 100);

        /**
         * @brief Create histogram of PSD values above energy threshold
         * @param psd_params Vector of PSD parameters
         * @param energy_threshold Energy threshold for histogram inclusion
         * @param n_bins Number of histogram bins
         * @return Generated histogram
         */
        static Histogram create_histogram(const std::vector<PSDParameters>& psd_params, double energy_threshold, int n_bins = 100);

        /**
         * @brief Set short integration gate duration
         * @param ns Short gate duration in nanoseconds
         */
        void set_short_gate(const double ns) { short_gate_ns_ = ns; }

        /**
         * @brief Set total integration gate duration
         * @param ns Total gate duration in nanoseconds
         */
        void set_total_gate(const double ns) { total_gate_ns_ = ns; }

        /**
         * @brief Set pre-trigger duration
         * @param ns Pre-trigger duration in nanoseconds
         */
        void set_pre_trigger(const double ns) { pre_trigger_ns_ = ns; }

        /**
         * @brief Get short integration gate duration
         * @return Short gate duration in nanoseconds
         */
        double get_short_gate() const { return short_gate_ns_; }

        /**
         * @brief Get total integration gate duration
         * @return Total gate duration in nanoseconds
         */
        double get_total_gate() const { return total_gate_ns_; }

        /**
         * @brief Get pre-trigger duration
         * @return Pre-trigger duration in nanoseconds
         */
        double get_pre_trigger() const { return pre_trigger_ns_; }

    private:
        double short_gate_ns_;
        double total_gate_ns_;
        double pre_trigger_ns_;

        /**
         * @brief Detect pile-up events in a waveform
         * @param waveform Waveform data
         * @param max_idx Index of maximum amplitude in the waveform
         * @param qtot Total integrated charge
         * @param psd_value Calculated PSD value
         * @return True if pile-up is detected, false otherwise
         */
        bool detect_pileup(const std::vector<double>& waveform, int max_idx, double qtot, double psd_value) const;

        /**
         * @brief Find peaks in histogram counts
         * @param counts Histogram counts
         * @param min_height Minimum height threshold for peak detection
         * @return Indices of detected peaks
         */
        static std::vector<int> find_peaks(const std::vector<int>& counts, double min_height = 0.15);

        /**
         * @brief Calculate Full Width at Half Maximum (FWHM) for a peak
         * @param counts Histogram counts
         * @param bin_centers Centers of histogram bins
         * @param peak_idx Index of the peak
         * @return Calculated FWHM value
         */
        static double calculate_fwhm(const std::vector<int>& counts, const std::vector<double>& bin_centers, int peak_idx);
    };
}