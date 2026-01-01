#pragma once

#include <vector>
#include <string>
#include <cstdint>

/// @brief Namespace for NE213 detector data analysis \namespace ne213
namespace ne213
{
    /**
     * @brief File type classification for measurement data \enum FileType
     *
     * Distinguishes between measurement files taken with the fusor running ('mit')
     * and background measurement files ('ohne'). This classification determines
     * whether PSD analysis should be performed.
     */
    enum class FileType
    {
        MIT,
        OHNE,
        UNKNOWN
    };

    /// @brief Teledyne LeCroy WaveAce 2024 sampling rate: 2 GS/s
    constexpr double SAMPLE_RATE_GS = 2.0;

    /// @brief Time per sample in nanoseconds (0.5 ns at 2 GS/s)
    constexpr double NS_PER_SAMPLE = 1.0 / SAMPLE_RATE_GS;

    /// @brief Samples per nanosecond (2.0 at 2 GS/s)
    constexpr double SAMPLES_PER_NS = 1.0 / NS_PER_SAMPLE;

    /// @brief Maximum number of samples to process per waveform
    constexpr size_t MAX_WAVEFORM_SAMPLES = 350;

    /**
     * @brief Convert time in nanoseconds to sample index
     * @param ns Time in nanoseconds
     * @return Corresponding sample index (integer)
     */
    inline int ns_to_sample(const double ns)
    {
        return static_cast<int>(ns * SAMPLES_PER_NS);
    }

    /**
     * @brief Convert sample index to time in nanoseconds
     * @param samples Sample index
     * @return Corresponding time in nanoseconds
     */
    inline double samples_to_ns(const int samples)
    {
        return samples / SAMPLES_PER_NS;
    }

    /**
     * @brief Convert sample index to time in nanoseconds (double overload)
     * @param samples Sample count as double
     * @return Corresponding time in nanoseconds
     */
    inline double samples_to_ns(const double samples)
    {
        return samples / SAMPLES_PER_NS;
    }

    /**
     * @brief Processed waveform data structure \struct WaveformData
     *
     * Contains baseline-corrected waveform data and extracted pulse characteristics
     * used for PSD analysis.
     */
    struct WaveformData
    {
        std::vector<double> data;
        double baseline{0.0};
        double max_amplitude{0.0};
        int max_index{0};
        double rise_time{0.0};
    };

    /**
     * @brief PSD analysis parameters for a single waveform \struct PSDParameters
     *
     * Contains the charge integration values and PSD parameter calculated using
     * Binda's Charge Comparison Method (Equation 8):
     *   PSD = (Qtot - Qshort) / Qtot
     *
     * Higher PSD values indicate neutrons (longer decay tail).
     * Lower PSD values indicate gamma rays (shorter decay tail).
     */
    struct PSDParameters
    {
        double psd_values{0.0};
        double qtot{0.0};
        double qshort{0.0};
        bool is_valid{false};
    };

    /**
     * @brief Figure of Merit (FOM) calculation result \struct FOMResult
     *
     * FOM quantifies neutron-gamma discrimination quality using Binda Equation (9):
     *   FOM = separation / (FWHM_gamma + FWHM_neutron)
     *
     * Reference values:
     *   - Binda (NE213 at JET): FOM = 2.25-2.35 (VERY GOOD)
     *   - Baselga (CCM baseline): FOM = 0.94
     *   - Baselga (best method): FOM = 1.04
     */
    struct FOMResult
    {
        double fom{0.0};
        double gamma_peak{0.0};
        double neutron_peak{0.0};
    };

    /**
     * @brief Histogram data structure for PSD distribution analysis \struct Histogram
     */
    struct Histogram
    {
        std::vector<double> bins;
        std::vector<int> counts;
        double min_value{0.0};
        double max_value{0.0};
    };

    /**
     * @brief Waveform file loader and classifier
     *
     * Handles loading waveform data from oscilloscope text files and classifies
     * files as measurement ('mit') or background ('ohne') based on filename.
     */
    class WaveformLoader
    {
    public:
        /**
         * @brief Load waveform data from a text file
         *
         * Reads waveforms from Teledyne LeCroy oscilloscope output format.
         * Handles CH1: prefix, 8-bit to signed conversion, and truncates
         * to MAX_WAVEFORM_SAMPLES (350 samples = 175 ns).
         *
         * @param filename Path to the input file
         * @return Vector of waveforms, each represented as a vector of doubles
         * @throws std::runtime_error if file cannot be opened
         */
        static std::vector<std::vector<double>> load_from_file(const std::string& filename);

        /**
         * @brief Classify file type based on filename
         *
         * Determines whether the file contains measurement data with the fusor
         * running ('mit') or background data ('ohne') by analyzing the filename.
         *
         * @param filename Path to the waveform data file
         * @return FileType::MIT if 'mit' found in filename,
         *         FileType::OHNE if 'ohne' found,
         *         FileType::UNKNOWN otherwise
         */
        static FileType classify_file_type(const std::string& filename);

        /**
         * @brief Convert FileType enum to human-readable string
         * @param type FileType enumeration value
         * @return String representation of the file type
         */
        static std::string file_type_to_string(FileType type);

    private:
        /**
         * @brief Parse a single line of waveform data
         *
         * Handles comma-separated values, CH1: prefix removal, and 8-bit
         * to signed conversion (values > 128 are converted to negative).
         *
         * @param line Input line of text
         * @return Parsed vector of doubles
         */
        static std::vector<double> parse_line(const std::string& line);
    };
}