#include "../include/WaveformProcessor.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <optional>
#include <iostream>
#include <iomanip>

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

double WaveformProcessor::calculate_baseline_std(const std::vector<double>& data, const size_t n_samples)
{
    if (data.empty()|| n_samples == 0) return 0.0;

    const size_t count = std::min(n_samples, data.size());
    double mean = 0.0;
    for (size_t i = 0; i < count; ++i)
    {
        mean += data[i];
    }
    mean /= static_cast<double>(count);

    double sum_eq = 0.0;
    for (size_t i = 0; i < count; ++i)
    {
        const double diff = data[i] - mean;
        sum_eq += diff * diff;
    }

    return std::sqrt(sum_eq / static_cast<double>(count));
}

BatchProcessingResult WaveformProcessor::process_batch(const std::vector<std::vector<double>>& raw_waveforms)
{
    BatchProcessingResult result;

    if (raw_waveforms.empty())
    {
        return result;
    }

    std::vector<std::vector<double>> inverted(raw_waveforms.size());
    for (size_t i = 0; i < raw_waveforms.size(); ++i)
    {
        inverted[i].resize(raw_waveforms[i].size());
        std::ranges::transform(raw_waveforms[i], inverted[i].begin(), [](const double x)
        {
           return -x;
        });
    }

    std::vector<double> baseline_stds;
    baseline_stds.reserve(inverted.size());
    for (const auto& wf : inverted)
    {
        baseline_stds.push_back(calculate_baseline_std(wf));
    }

    std::vector<double> sorted_stds = baseline_stds;
    std::ranges::sort(sorted_stds);
    const double median_std = sorted_stds[sorted_stds.size() / 2];
    const double noise_threshold = median_std * 3.0;

    std::cout << "Baseline noise analysis:\n";
    std::cout << "  Median baseline std: " << std::fixed << std::setprecision(2) << median_std << "\n";
    std::cout << "  Adaptive noise threshold: " << noise_threshold << "\n";

    for (size_t i = 0; i  < inverted.size(); ++i)
    {
        const double bl_std = baseline_stds[i];

        if (bl_std > noise_threshold)
        {
            result.rejected.push_back({i, "high_baseline_noise", bl_std, noise_threshold});
            result.noise_rejections++;
            continue;
        }

        const double baseline = calculate_baseline(inverted[i]);
        std::vector<double> corrected = apply_baseline_correction(inverted[i], baseline);

        const size_t check_end = std::min(corrected.size(), static_cast<size_t>(100));
        const double min_pretrigger = *std::min_element(corrected.begin(), corrected.begin() + static_cast<long>(check_end));
        const double spike_threshold = -(6.0 * bl_std + 100.0);

        if (min_pretrigger < spike_threshold)
        {
            result.rejected.push_back({i, "negative_spike", min_pretrigger, spike_threshold});
            result.spike_rejections++;
            continue;
        }

        WaveformData waveform;
        waveform.data = std::move(corrected);
        waveform.baseline = baseline;

        auto [max_amp, max_idx] = find_peak(waveform.data);
        waveform.max_amplitude = max_amp;
        waveform.max_index = max_idx;
        waveform.rise_time = calculate_rise_time(waveform.data, max_amp);

        result.waveforms.push_back(std::move(waveform));
    }

    const size_t total = raw_waveforms.size();
    const size_t rejected = result.rejected.size();
    std::cout << "  Waveforms rejected: " << rejected << " / "
                << total << " (" << std::fixed << std::setprecision(1) << (100.0 * static_cast<double>(rejected) / static_cast<double>(total)) << "%)\n";

    if (!result.rejected.empty())
    {
        std::cout << "      - High baseline noise: " << result.noise_rejections << "\n";
        std::cout << "      - Negative spikes: " << result.spike_rejections << "\n";
    }

    return result;
}
