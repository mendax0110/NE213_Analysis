#pragma once

#include "WaveformData.h"
#include <vector>
#include <string>
#include <sstream>

/// @brief Namespace for NE213 detector data analysis \namespace ne213
namespace ne213
{
    /**
     * @brief Complete statistical report for NE213 analysis \struct StatisticsReport
     *
     * Contains all computed statistics from a measurement session including
     * event counts, pulse characteristics, PSD statistics, and FOM results.
     */
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

        FileType file_type{FileType::UNKNOWN};
    };

    /**
     * @brief Statistical analysis and report generation
     *
     * Static class providing comprehensive statistical analysis of NE213
     * measurement data and formatted report output.
     */
    class Statistics
    {
    public:
        /**
         * @brief Generate comprehensive statistics report
         *
         * Computes all statistics from processed waveforms and PSD parameters
         * including event counts, pulse characteristics, PSD distribution
         * statistics, and Figure of Merit.
         *
         * @param waveforms Vector of processed waveform data
         * @param psd_params Vector of corresponding PSD parameters
         * @param short_gate_ns Short gate duration used for PSD calculation
         * @param total_gate_ns Total gate duration used for PSD calculation
         * @param energy_threshold Minimum Qtot for FOM calculation
         * @param file_type Classification of the input file (MIT/OHNE/UNKNOWN)
         * @return Complete statistics report
         */
        static StatisticsReport generate_report(
            const std::vector<WaveformData>& waveforms,
            const std::vector<PSDParameters>& psd_params,
            double short_gate_ns,
            double total_gate_ns,
            double energy_threshold = 10.0,
            FileType file_type = FileType::MIT
        );

        /**
         * @brief Print formatted statistics report to console
         *
         * Outputs a formatted report to stdout with sections for detector
         * configuration, pulse characteristics, Binda CCM analysis, and
         * discrimination performance assessment.
         *
         * @param report Statistics report to print
         */
        static void print_report(const StatisticsReport& report);

        /**
         * @brief Generate report as a string
         *
         * Creates the same formatted report as print_report but returns
         * it as a string for display in GUI or file output.
         *
         * @param report Statistics report to format
         * @return Formatted report as string
         */
        static std::string generate_report_string(const StatisticsReport& report);

        /**
         * @brief Get FOM quality assessment text
         *
         * Returns a text description of the FOM quality based on reference values:
         * - FOM > 2.0: "VERY GOOD (matching Binda JET results)"
         * - FOM > 1.0: "GOOD discrimination capability"
         * - FOM > 0.7: "MODERATE discrimination"
         * - Otherwise: "POOR discrimination"
         *
         * @param fom Figure of Merit value
         * @return Quality assessment string
         */
        static std::string get_fom_quality_text(double fom);

    private:
        /**
         * @brief Calculate arithmetic mean of values
         * @param values Input vector of values
         * @return Mean value, or 0 if empty
         */
        static double calculate_mean(const std::vector<double>& values);

        /**
         * @brief Calculate sample standard deviation
         * @param values Input vector of values
         * @param mean Pre-calculated mean value
         * @return Standard deviation, or 0 if fewer than 2 values
         */
        static double calculate_std(const std::vector<double>& values, double mean);
    };
}