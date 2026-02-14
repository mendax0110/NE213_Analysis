#pragma once

#include "WaveformData.h"
#include <vector>
#include <optional>

/// @brief Namespace for NE213 detector data analysis \namespace ne213
namespace ne213
{
    /**
     * @brief Waveform signal processor for NE213 scintillator pulses
     *
     * Static class providing signal processing functions for preparing raw
     * oscilloscope waveforms for PSD analysis. Handles inversion, baseline
     * subtraction, and pulse characteristic extraction.
     */
    class WaveformProcessor
    {
    public:
        /**
         * @brief Process a raw waveform for PSD analysis
         *
         * Performs the complete processing pipeline:
         * 1. Invert waveform (NE213 pulses are negative-going)
         * 2. Calculate and subtract baseline (mean of first 50 samples)
         * 3. Find peak amplitude and index
         * 4. Calculate rise time (10% to 90% of peak)
         *
         * @param raw_data Vector of raw waveform samples from oscilloscope
         * @return Processed waveform data with extracted characteristics
         */
        static WaveformData process_waveform(const std::vector<double>& raw_data);

        /**
         * @brief Process a batch of raw waveforms
         *
         * Applies the processing pipeline to each waveform in the input vector.
         * Collects results and counts of rejected events due to noise or spikes.
         *
         * @param raw_waveforms Vector of raw waveforms (each a vector of samples)
         * @return Batch processing result containing processed waveforms and rejection info
         */
        static BatchProcessingResult process_batch(const std::vector<std::vector<double>>& raw_waveforms);

        /**
         * @brief Calculate baseline from waveform pre-trigger region
         *
         * Computes the mean of the first n_samples as the baseline level.
         * For NE213 signals, the pre-trigger region should be stable.
         *
         * @param data Vector of waveform samples
         * @param n_samples Number of initial samples to use (default: 50)
         * @return Calculated baseline value
         */
        static double calculate_baseline(const std::vector<double>& data, size_t n_samples = 50);

        /**
         * @brief Apply baseline correction
         *
         * Subtracts baseline from all samples. Does not clip negative values,
         * matching the Python implementation behavior.
         *
         * @param data Vector of waveform samples
         * @param baseline Baseline value to subtract
         * @return Baseline-corrected waveform data
         */
        static std::vector<double> apply_baseline_correction(const std::vector<double>& data, double baseline);

        /**
         * @brief Calculate pulse rise time (10% to 90% of peak)
         *
         * Measures the time for the pulse to rise from 10% to 90% of its
         * peak amplitude. This is a characteristic of the scintillation
         * process and can help identify pulse quality.
         *
         * @param data Vector of baseline-corrected waveform samples
         * @param max_amplitude Maximum amplitude of the waveform
         * @return Rise time in sample units (multiply by NS_PER_SAMPLE for ns)
         */
        static double calculate_rise_time(const std::vector<double>& data, double max_amplitude);

        /**
         * @brief Find peak amplitude and its index in the waveform
         *
         * Locates the maximum value in the waveform and returns both
         * the amplitude and sample index where it occurs.
         *
         * @param data Vector of waveform samples
         * @return Pair of (peak amplitude, sample index)
         */
        static std::pair<double, int> find_peak(const std::vector<double>& data);

    private:
        /**
         * @brief Find first crossing of a threshold level
         *
         * Scans the waveform from the beginning to find the first sample
         * where the amplitude exceeds the specified threshold.
         *
         * @param data Vector of waveform samples
         * @param threshold Threshold value to find crossing
         * @return Optional index of threshold crossing, or std::nullopt if not found
         */
        static std::optional<size_t> find_threshold_crossing(const std::vector<double>& data, double threshold);

        /**
         * @brief Calculate standard devication of baseline region
         * @param data Vector of waveform samples
         * @param n_samples Number of initial samples to use for baseline std calculation (default: 50)
         * @return Population standard deviation of the baseline region.
         */
        static double calculate_baseline_std(const std::vector<double>& data, size_t n_samples = 50);
    };
}