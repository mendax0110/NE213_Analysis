#pragma once

#include "WaveformData.h"
#include <vector>
#include <optional>

/// @brief Namespace for NE213 detector data analysis \namespace ne213
namespace ne213
{
    /// @brief Class for processing waveforms \class WaveformProcessor
    class WaveformProcessor
    {
    public:
        /**
         * @brief Process raw waveform data
         * @param raw_data Vector of raw waveform samples
         * @return Processed waveform data
         */
        static WaveformData process_waveform(const std::vector<double>& raw_data);

        /**
         * @brief Calculate baseline from waveform data
         * @param data Vector of waveform samples
         * @param n_samples Number of initial samples to use for baseline calculation
         * @return Calculated baseline value
         */
        static double calculate_baseline(const std::vector<double>& data, size_t n_samples = 50);

        /**
         * @brief Apply baseline correction to waveform data
         * @param data Vector of waveform samples
         * @param baseline Baseline value to subtract
         * @return Baseline-corrected waveform data
         */
        static std::vector<double> apply_baseline_correction(const std::vector<double>& data, double baseline);

        /**
         * @brief Calculate rise time of the waveform
         * @param data Vector of waveform samples
         * @param max_amplitude Maximum amplitude of the waveform
         * @return Calculated rise time in sample units
         */
        static double calculate_rise_time(const std::vector<double>& data, double max_amplitude);

        /**
         * @brief Find the peak amplitude and its index in the waveform
         * @param data Vector of waveform samples
         * @return Pair of peak amplitude and its index
         */
        static std::pair<double, int> find_peak(const std::vector<double>& data);

    private:
        /**
         * @brief Find the first crossing of a specified threshold in the waveform
         * @param data Vector of waveform samples
         * @param threshold Threshold value to find crossing
         * @return Optional index of the threshold crossing, or std::nullopt if not found
         */
        static std::optional<size_t> find_threshold_crossing(const std::vector<double>& data, double threshold);
    };
}