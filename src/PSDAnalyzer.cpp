#include "../include/PSDAnalyzer.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>

using namespace ne213;

PSDAnalyzer::PSDAnalyzer(const double short_gate_ns, const double total_gate_ns, const double pre_trigger_ns)
    : short_gate_ns_(short_gate_ns)
    , total_gate_ns_(total_gate_ns)
    , pre_trigger_ns_(pre_trigger_ns)
{

}

PSDParameters PSDAnalyzer::calculate_parameters(const WaveformData& waveform) const
{
    PSDParameters params;

    if (waveform.max_index < 0 || waveform.max_index >= static_cast<int>(waveform.data.size()) - 10)
    {
        return params;
    }

    const int short_gate_samples = ns_to_sample(short_gate_ns_);
    const int total_gate_samples = ns_to_sample(total_gate_ns_);
    const int pre_trigger_samples = ns_to_sample(pre_trigger_ns_);

    const int start_idx = std::max(0, waveform.max_index - pre_trigger_samples);

    const int total_end = std::min(
        static_cast<int>(waveform.data.size()),
        start_idx + total_gate_samples
    );
    params.qtot = std::accumulate(
        waveform.data.begin() + start_idx,
        waveform.data.begin() + total_end,
        0.0
    );

    const int short_end = std::min(
        static_cast<int>(waveform.data.size()),
        start_idx + short_gate_samples
    );
    params.qshort = std::accumulate(
        waveform.data.begin() + start_idx,
        waveform.data.begin() + short_end,
        0.0
    );

    if (params.qtot > 0)
    {
        params.psd_values = (params.qtot - params.qshort) / params.qtot;
    }

    params.is_valid = !detect_pileup(
        waveform.data,
        waveform.max_index,
        params.qtot,
        params.psd_values
    );

    return params;
}

bool PSDAnalyzer::detect_pileup(const std::vector<double>& waveform, const int max_idx, double qtot, const double psd_value) const
{
    if (max_idx < 20 || max_idx >= static_cast<int>(waveform.size()) - 10)
    {
        return true;
    }

    if (psd_value > 0.5)
    {
        return true;
    }

    const int rising_start = std::max(0, max_idx - 20);
    std::vector<double> rising_edge(
        waveform.begin() + rising_start,
        waveform.begin() + max_idx
    );

    if (rising_edge.size() > 10)
    {
        const double max_rise = *std::max_element(rising_edge.begin(), rising_edge.end());
        const double threshold = 0.25 * max_rise;

        int peak_count = 0;
        for (size_t i = 1; i < rising_edge.size() - 1; ++i)
        {
            if (rising_edge[i] >= threshold && rising_edge[i] > rising_edge[i-1] && rising_edge[i] > rising_edge[i+1])
            {
                peak_count++;
            }
        }

        if (peak_count > 1)
        {
            return true;
        }
    }

    if (max_idx + 30 < static_cast<int>(waveform.size()))
    {
        for (int i = max_idx + 1; i < max_idx + 30; ++i)
        {
            if (waveform[i] > waveform[i-1] * 1.15)
            {
                return true;
            }
        }
    }

    return false;
}

FOMResult PSDAnalyzer::calculate_fom(const std::vector<PSDParameters>& psd_params, const double energy_threshold)
{
    FOMResult result;

    std::vector<double> valid_psd;
    for (const auto& param : psd_params)
    {
        if (param.is_valid && param.qtot > energy_threshold)
        {
            valid_psd.push_back(param.psd_values);
        }
    }

    if (valid_psd.size() < 20)
    {
        return result;
    }

    const auto hist = create_histogram(psd_params, energy_threshold, 100);
    auto peaks = find_peaks(hist.counts);

    if (peaks.size() < 2)
    {
        return result;
    }

    std::ranges::sort(peaks);
    const int gamma_peak_idx = peaks[0];
    const int neutron_peak_idx = peaks[peaks.size() - 1];

    result.gamma_peak = hist.bins[gamma_peak_idx];
    result.neutron_peak = hist.bins[neutron_peak_idx];

    const double separation = std::abs(result.neutron_peak - result.gamma_peak);
    const double gamma_fwhm = calculate_fwhm(hist.counts, hist.bins, gamma_peak_idx);
    const double neutron_fwhm = calculate_fwhm(hist.counts, hist.bins, neutron_peak_idx);

    result.fom = separation / (gamma_fwhm + neutron_fwhm);

    return result;
}

Histogram PSDAnalyzer::create_histogram(const std::vector<PSDParameters>& psd_params, const double energy_threshold, const int n_bins)
{
    Histogram hist;

    std::vector<double> valid_psd;
    for (const auto& param : psd_params)
    {
        if (param.is_valid && param.qtot > energy_threshold)
        {
            valid_psd.push_back(param.psd_values);
        }
    }

    if (valid_psd.empty())
    {
        return hist;
    }

    hist.min_value = *std::ranges::min_element(valid_psd);
    hist.max_value = *std::ranges::max_element(valid_psd);

    const double bin_width = (hist.max_value - hist.min_value) / n_bins;
    hist.bins.resize(n_bins);
    hist.counts.resize(n_bins, 0);

    for (int i = 0; i < n_bins; ++i)
    {
        hist.bins[i] = hist.min_value + (i + 0.5) * bin_width;
    }

    for (const double psd : valid_psd)
    {
        int bin = static_cast<int>((psd - hist.min_value) / bin_width);
        bin = std::clamp(bin, 0, n_bins - 1);
        hist.counts[bin]++;
    }

    return hist;
}

std::vector<int> PSDAnalyzer::find_peaks(const std::vector<int>& counts, double min_height)
{
    std::vector<int> peaks;
    if (counts.size() < 3) return peaks;

    const int max_count = *std::ranges::max_element(counts);
    const int height_threshold = static_cast<int>(max_count * min_height);
    const int prominence_threshold = static_cast<int>(max_count * 0.1);
    constexpr int min_distance = 10;

    std::vector<std::pair<int, int>> candidates;
    for (size_t i = 1; i < counts.size() - 1; ++i)
    {
        if (counts[i] > height_threshold && counts[i] > counts[i-1] && counts[i] > counts[i+1])
        {
            candidates.push_back({static_cast<int>(i), counts[i]});
        }
    }

    std::ranges::sort(candidates,[](const auto& a, const auto& b)
    {
        return a.second > b.second;
    });

    for (const auto& [idx, height] : candidates)
    {
        int left_min = height;
        for (int j = idx - 1; j >= 0 && j >= idx - 20; --j)
        {
            left_min = std::min(left_min, counts[j]);
        }

        int right_min = height;
        for (int j = idx + 1; j < static_cast<int>(counts.size()) && j <= idx + 20; ++j)
        {
            right_min = std::min(right_min, counts[j]);
        }

        const int prominence = height - std::max(left_min, right_min);
        if (prominence < prominence_threshold) continue;

        bool too_close = false;
        for (const int existing_peak : peaks)
        {
            if (std::abs(idx - existing_peak) < min_distance)
            {
                too_close = true;
                break;
            }
        }

        if (!too_close)
        {
            peaks.push_back(idx);
        }
    }

    std::ranges::sort(peaks);

    return peaks;
}

double PSDAnalyzer::calculate_fwhm(const std::vector<int>& counts, const std::vector<double>& bin_centers, const int peak_idx)
{
    if (peak_idx < 0 || peak_idx >= static_cast<int>(counts.size()))
    {
        return 0.01;
    }

    const double peak_height = counts[peak_idx];
    const double half_max = peak_height / 2.0;

    int left_idx = peak_idx;
    while (left_idx > 0 && counts[left_idx] > half_max)
    {
        left_idx--;
    }

    int right_idx = peak_idx;
    while (right_idx < static_cast<int>(counts.size()) - 1 && counts[right_idx] > half_max)
    {
        right_idx++;
    }

    if (left_idx < peak_idx && peak_idx < right_idx)
    {
        const double fwhm = bin_centers[right_idx] - bin_centers[left_idx];
        return std::max(fwhm, 0.001);
    }

    return 0.01;
}