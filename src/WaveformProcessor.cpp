#include "../include/WaveformProcessor.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <optional>

using namespace ne213;

WaveformData WaveformProcessor::process_waveform(const std::vector<double>& raw_data)
{
    WaveformData waveform;

    std::vector<double> inverted(raw_data.size());
    std::ranges::transform(raw_data, inverted.begin(),[](const double x)
    {
        return -x;
    });

    waveform.baseline = calculate_baseline(inverted);
    waveform.data = apply_baseline_correction(inverted, waveform.baseline);

    auto [max_amp, max_idx] = find_peak(waveform.data);
    waveform.max_amplitude = max_amp;
    waveform.max_index = max_idx;
    waveform.rise_time = calculate_rise_time(waveform.data, max_amp);

    return waveform;
}

double WaveformProcessor::calculate_baseline(const std::vector<double>& data, const size_t n_samples)
{
    if (data.empty() || n_samples == 0) return 0.0;

    const size_t count = std::min(n_samples, data.size());
    const double sum = std::accumulate(data.begin(), data.begin() + count, 0.0);
    return sum / static_cast<double>(count);
}

std::vector<double> WaveformProcessor::apply_baseline_correction(const std::vector<double>& data, double baseline)
{
    std::vector<double> corrected(data.size());
    std::ranges::transform(data, corrected.begin(),[baseline](double x)
    {
        return std::max(0.0, x - baseline);
    });
    return corrected;
}

double WaveformProcessor::calculate_rise_time(const std::vector<double>& data, const double max_amplitude)
{
    if (max_amplitude <= 0) return 0.0;

    const double threshold_10 = 0.1 * max_amplitude;
    const double threshold_90 = 0.9 * max_amplitude;

    const auto idx_10 = find_threshold_crossing(data, threshold_10);
    const auto idx_90 = find_threshold_crossing(data, threshold_90);

    if (idx_10.has_value() && idx_90.has_value() && *idx_90 > *idx_10)
    {
        return static_cast<double>(*idx_90 - *idx_10);
    }
    return 0.0;
}

std::pair<double, int> WaveformProcessor::find_peak(const std::vector<double>& data)
{
    if (data.empty()) return {0.0, 0};

    auto max_it = std::ranges::max_element(data);
    return {*max_it, static_cast<int>(std::distance(data.begin(), max_it))};
}

std::optional<size_t> WaveformProcessor::find_threshold_crossing(const std::vector<double>& data, const double threshold)
{
    for (size_t i = 0; i < data.size(); ++i)
    {
        if (data[i] > threshold)
        {
            return i;
        }
    }
    return std::nullopt;
}