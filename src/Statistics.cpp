#include "../include/Statistics.h"
#include "../include/PSDAnalyzer.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <sstream>

using namespace ne213;

StatisticsReport Statistics::generate_report(
    const std::vector<WaveformData>& waveforms,
    const std::vector<PSDParameters>& psd_params,
    const double short_gate_ns,
    const double total_gate_ns,
    const double energy_threshold,
    const FileType file_type)
{
    StatisticsReport report;

    report.total_events = waveforms.size();
    report.short_gate_ns = short_gate_ns;
    report.total_gate_ns = total_gate_ns;
    report.file_type = file_type;

    std::vector<double> valid_psd;
    std::vector<double> valid_qtot;
    std::vector<double> valid_qshort;
    std::vector<double> amplitudes;
    std::vector<double> rise_times;

    for (size_t i = 0; i < waveforms.size(); ++i)
    {
        amplitudes.push_back(waveforms[i].max_amplitude);
        rise_times.push_back(waveforms[i].rise_time * NS_PER_SAMPLE);

        if (i < psd_params.size() && psd_params[i].is_valid)
        {
            report.valid_events++;
            valid_psd.push_back(psd_params[i].psd_values);
            valid_qtot.push_back(psd_params[i].qtot);
            valid_qshort.push_back(psd_params[i].qshort);
        }
    }

    report.pileup_events = report.total_events - report.valid_events;
    report.mean_amplitude = calculate_mean(amplitudes);
    report.mean_rise_time_ns = calculate_mean(rise_times);

    if (!valid_psd.empty())
    {
        report.psd_min = *std::ranges::min_element(valid_psd);
        report.psd_max = *std::ranges::max_element(valid_psd);
        report.psd_mean = calculate_mean(valid_psd);
        report.psd_std = calculate_std(valid_psd, report.psd_mean);

        std::vector<double> sorted_psd = valid_psd;
        std::ranges::sort(sorted_psd);
        if (sorted_psd.size() % 2 == 0)
        {
            report.psd_median = (sorted_psd[sorted_psd.size()/2 - 1] + sorted_psd[sorted_psd.size()/2]) / 2.0;
        }
        else
        {
            report.psd_median = sorted_psd[sorted_psd.size()/2];
        }

        report.qtot_mean = calculate_mean(valid_qtot);
        report.qshort_mean = calculate_mean(valid_qshort);

        if (report.qtot_mean > 0)
        {
            report.slow_component_fraction = (report.qtot_mean - report.qshort_mean) / report.qtot_mean;
        }
    }

    if (file_type == FileType::MIT)
    {
        report.fom_result = PSDAnalyzer::calculate_fom(psd_params, energy_threshold);

        if (!valid_psd.empty() && report.valid_events > 0)
        {
            for (const auto& param : psd_params)
            {
                if (param.is_valid && param.psd_values > report.psd_median)
                {
                    report.estimated_neutron_count++;
                }
            }
            report.neutron_fraction = static_cast<double>(report.estimated_neutron_count) / static_cast<double>(report.valid_events);
        }
    }

    return report;
}

std::string Statistics::get_fom_quality_text(const double fom)
{
    if (fom > 2.0)
    {
        return "VERY GOOD (matching Binda JET results)";
    }
    if (fom > 1.0)
    {
        return "GOOD discrimination capability";
    }
    if (fom > 0.7)
    {
        return "MODERATE discrimination";
    }
    return "POOR discrimination";
}

std::string Statistics::generate_report_string(const StatisticsReport& report)
{
    std::ostringstream ss;

    ss << std::string(70, '=') << "\n";
    ss << "NE213 NEUTRON DETECTION ANALYSIS - Farnsworth Fusion Reactor\n";
    ss << "File type: " << WaveformLoader::file_type_to_string(report.file_type) << "\n";
    ss << "Binda Charge Comparison Method Implementation\n";
    ss << std::string(70, '=') << "\n";
    ss << "Total events captured: " << report.total_events << "\n";

    if (report.file_type == FileType::MIT)
    {
        const double valid_pct = report.total_events > 0 ? (100.0 * report.valid_events / report.total_events) : 0.0;
        const double pileup_pct = report.total_events > 0 ? (100.0 * report.pileup_events / report.total_events) : 0.0;

        ss << std::fixed << std::setprecision(1);
        ss << "Valid events (pile-up rejected): " << report.valid_events << " (" << valid_pct << "%)\n";
        ss << "Pile-up events rejected: " << report.pileup_events << " (" << pileup_pct << "%)\n\n";
    }
    else
    {
        ss << "Background measurement - no pile-up rejection applied\n\n";
    }

    ss << "Detector Configuration:\n";
    ss << "  Type: NE213 Liquid Scintillator\n";
    ss << "  Target reaction: D-D fusion\n";
    ss << "  Sampling rate: " << SAMPLE_RATE_GS << " GS/s (" << NS_PER_SAMPLE << " ns/sample)\n\n";

    ss << "Pulse Characteristics:\n";
    ss << std::setprecision(2);
    ss << "  Average amplitude: " << report.mean_amplitude << "\n";
    ss << "  Average rise time: " << report.mean_rise_time_ns << " ns\n\n";

    if (report.file_type == FileType::MIT)
    {
        ss << "Binda Charge Comparison Method (Equation 8):\n";
        ss << "   PSD = (Qtot - Qshort) / Qtot\n";
        ss << "  Short gate: " << report.short_gate_ns << " ns\n";
        ss << "  Total gate: " << report.total_gate_ns << " ns\n";
        ss << std::setprecision(4);
        ss << "  PSD range: [" << report.psd_min << ", " << report.psd_max << "]\n";
        ss << "  PSD mean: " << report.psd_mean << "\n";
        ss << "  PSD std: " << report.psd_std << "\n\n";

        ss << "Integrated Charges:\n";
        ss << std::setprecision(1);
        ss << "  Qtot mean: " << report.qtot_mean << "\n";
        ss << "  Qshort mean: " << report.qshort_mean << "\n";
        ss << "  Slow component fraction: " << (report.slow_component_fraction * 100.0) << "%\n\n";

        ss << "Discrimination Performance (Binda Section 6.2):\n";
        ss << std::setprecision(4);
        ss << "  Figure of Merit (FOM): " << report.fom_result.fom << "\n";
        ss << "  Reference values:\n";
        ss << "    - Binda (NE213 at JET): FOM = 2.25-2.35 (VERY GOOD)\n";
        ss << "    - Baselga (CCM baseline): FOM = 0.94\n";
        ss << "    - Baselga (best method): FOM = 1.04\n";
        ss << "  [RESULT] " << get_fom_quality_text(report.fom_result.fom) << "\n";

        if (report.fom_result.fom > 0.7)
        {
            ss << "  Estimated neutron events: " << report.estimated_neutron_count
               << " (" << std::setprecision(1) << (report.neutron_fraction * 100.0) << "%)\n";
        }

        if (report.fom_result.gamma_peak > 0 && report.fom_result.neutron_peak > 0)
        {
            ss << "\nPeak Positions:\n";
            ss << std::setprecision(4);
            ss << "  Gamma peak (Y): " << report.fom_result.gamma_peak << "\n";
            ss << "  Neutron peak (n): " << report.fom_result.neutron_peak << "\n";
            ss << "  Separation: " << (report.fom_result.neutron_peak - report.fom_result.gamma_peak) << "\n";
            ss << std::setprecision(2);
            ss << "  Peak ratio: " << (report.fom_result.neutron_peak / report.fom_result.gamma_peak) << "\n";
        }
    }
    else
    {
        ss << "Background Measurement Summary:\n";
        ss << "  This data represents background noise (cosmic rays, natural radioactivity, etc.)\n";
        ss << "  Use this for noise characterization and background subtraction\n";
        ss << "  No neutron/gamma discrimination analysis performed on background data\n";
    }

    ss << "\n";
    return ss.str();
}

void Statistics::print_report(const StatisticsReport& report)
{
    std::cout << generate_report_string(report);
}

double Statistics::calculate_mean(const std::vector<double>& values)
{
    if (values.empty()) return 0.0;
    return std::accumulate(values.begin(), values.end(), 0.0) / static_cast<double>(values.size());
}

double Statistics::calculate_std(const std::vector<double>& values, const double mean)
{
    if (values.size() <= 1) return 0.0;

    double sum_sq = 0.0;
    for (const double val : values)
    {
        const double diff = val - mean;
        sum_sq += diff * diff;
    }
    return std::sqrt(sum_sq / static_cast<double>(values.size() - 1));
}
