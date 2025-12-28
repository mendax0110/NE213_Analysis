#pragma once

#include "WaveformData.h"
#include <vector>
#include <string>

/// @brief Namespace for NE213 detector data analysis \namespace ne213
namespace ne213
{
    /// @brief Structure to hold statistical report data \struct StatisticsReport
    struct StatisticsReport
    {
        size_t total_events{0};
        size_t valid_events{0};
        size_t pileup_events{0};

        double mean_amplitude{0.0};
        double mean_rise_time_ns{0.0};

        double psd_min{0.0};
        double psd_max{0.0};
        double psd_mean{0.0};
        double psd_std{0.0};
        double psd_median{0.0};

        double qtot_mean{0.0};
        double qshort_mean{0.0};
        double slow_component_fraction{0.0};

        FOMResult fom_result;
        size_t estimated_neutron_count{0};
        double neutron_fraction{0.0};

        double short_gate_ns{0.0};
        double total_gate_ns{0.0};
    };

    /// @brief Class for generating and printing statistical reports \class Statistics
    class Statistics
    {
    public:
        /**
         * @brief Generate comprehensive statistics report
         * @param waveforms Processed waveforms
         * @param psd_params PSD parameters
         * @param short_gate_ns Short integration gate duration
         * @param total_gate_ns Total integration gate duration
         * @param energy_threshold Energy threshold for FOM calculation
         * @return Statistics report
         */
        static StatisticsReport generate_report(
            const std::vector<WaveformData>& waveforms,
            const std::vector<PSDParameters>& psd_params,
            double short_gate_ns,
            double total_gate_ns,
            double energy_threshold = 10.0
        );

        /**
         * @brief Print formatted statistics report to console
         * @param report Statistics report
         */
        static void print_report(const StatisticsReport& report);

    private:
        /**
         * @brief Calculate mean of a vector of values
         * @param values Input values
         * @return Mean value
         */
        static double calculate_mean(const std::vector<double>& values);

        /**
         * @brief Calculate standard deviation of a vector of values
         * @param values Input values
         * @param mean Pre-calculated mean value
         * @return Standard deviation
         */
        static double calculate_std(const std::vector<double>& values, double mean);
    };

} // namespace ne213