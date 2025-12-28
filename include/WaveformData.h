#pragma once

#include <vector>
#include <string>
#include <cstdint>

/// @brief Namespace for NE213 detector data analysis \namespace ne213
namespace ne213
{
    /// @brief Sample rate and conversion constants
    constexpr double SAMPLE_RATE_GS = 2.0;
    constexpr double NS_PRE_SAMPLE = 1.0 / SAMPLE_RATE_GS;
    constexpr double SAMPLES_PER_NS = 1.0 / NS_PRE_SAMPLE;

    /**
     * @brief Convert nanoseconds to sample index
     * @param ns Time in nanoseconds
     * @return Corresponding sample index
     */
    inline int ns_to_sample(const double ns)
    {
        return static_cast<int>(ns * SAMPLES_PER_NS);
    }

    /**
     * @brief Convert sample index to nanoseconds
     * @param samples Sample index
     * @return Corresponding time in nanoseconds
     */
    inline double samples_to_ns(const int samples)
    {
        return samples / SAMPLES_PER_NS;
    }

    /// @brief Structure to hold processed waveform data \struct WaveformData
    struct WaveformData
    {
        std::vector<double> data;
        double baseline{0.0};
        double max_amplitude{0.0};
        int max_index{0};
        double rise_time{0.0};
    };

    /// @brief Structure to hold PSD parameters \struct PSDParameters
    struct PSDParameters
    {
        double psd_values{0.0};
        double qtot{0.0};
        double qshort{0.0};
        bool is_valid{false};
    };

    /// @brief Structure to hold Figure of Merit (FOM) results \struct FOMResult
    struct FOMResult
    {
        double fom{0.0};
        double gamma_peak{0.0};
        double neutron_peak{0.0};
    };

    /// @brief Structure to hold histogram data \struct Histogram
    struct Histogram
    {
        std::vector<double> bins;
        std::vector<int> counts;
        double min_value{0.0};
        double max_value{0.0};
    };

    /// @brief Class for loading waveform data from files \class WaveformLoader
    class WaveformLoader
    {
    public:
        /**
         * @brief Load waveform data from a text file
         * @param filename Path to the input file
         * @return Vector of waveforms, each represented as a vector of doubles
         */
        static std::vector<std::vector<double>> load_from_file(const std::string& filename);

    private:
        /**
         * @brief Parse a line of text into a vector of doubles
         * @param line Input line of text
         * @return Parsed vector of doubles
         */
        static std::vector<double> parse_line(const std::string& line);
    };
}