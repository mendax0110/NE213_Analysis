#include "../include/Statistics.h"
#include "../include/PSDAnalyzer.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <numeric>
#include <algorithm>


using namespace ne213;


StatisticsReport Statistics::generate_report(
    const std::vector<WaveformData>& waveforms,
    const std::vector<PSDParameters>& psd_params,
    const double short_gate_ns,
    const double total_gate_ns,
    const double energy_threshold)
{
    StatisticsReport report;

    report.total_events = waveforms.size();
    report.short_gate_ns = short_gate_ns;
    report.total_gate_ns = total_gate_ns;

    std::vector<double> valid_psd;
    std::vector<double> valid_qtot;
    std::vector<double> valid_qshort;
    std::vector<double> amplitudes;
    std::vector<double> rise_times;

    for (size_t i = 0; i < waveforms.size(); ++i)
    {
        amplitudes.push_back(waveforms[i].max_amplitude);
        rise_times.push_back(waveforms[i].rise_time * NS_PRE_SAMPLE);

        if (psd_params[i].is_valid)
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

    report.fom_result = PSDAnalyzer::calculate_fom(psd_params, energy_threshold);

    if (!valid_psd.empty())
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

    return report;
}

void Statistics::print_report(const StatisticsReport& report)
{
    std::cout << std::string(70, '=') << "\n";
    std::cout << "NE213 NEUTRON DETECTION ANALYSIS - Farnsworth Fusion Reactor\n";
    std::cout << "Binda Charge Comparison Method Implementation\n";
    std::cout << std::string(70, '=') << "\n";
    std::cout << "Total events captured: " << report.total_events << "\n";
    std::cout << "Valid events (pile-up rejected): " << report.valid_events
              << " (" << std::fixed << std::setprecision(1)
              << (100.0 * report.valid_events / report.total_events) << "%)\n";
    std::cout << "Pile-up events rejected: " << report.pileup_events
              << " (" << (100.0 * report.pileup_events / report.total_events) << "%)\n\n";

    std::cout << "Detector Configuration:\n";
    std::cout << "  Type: NE213 Liquid Scintillator\n";
    std::cout << "  Target reaction: D-D fusion\n";
    std::cout << "  Sampling rate: " << SAMPLE_RATE_GS << " GS/s ("
              << NS_PRE_SAMPLE << " ns/sample)\n\n";

    std::cout << "Pulse Characteristics:\n";
    std::cout << "  Average amplitude: " << std::setprecision(2)
              << report.mean_amplitude << "\n";
    std::cout << "  Average rise time: " << report.mean_rise_time_ns << " ns\n\n";

    std::cout << "Binda Charge Comparison Method (Equation 8):\n";
    std::cout << "   PSD = (Qtot - Qshort) / Qtot\n";
    std::cout << "  Short gate: " << report.short_gate_ns << " ns\n";
    std::cout << "  Total gate: " << report.total_gate_ns << " ns\n";
    std::cout << "  PSD range: [" << std::setprecision(4) << report.psd_min
              << ", " << report.psd_max << "]\n";
    std::cout << "  PSD mean: " << report.psd_mean << "\n";
    std::cout << "  PSD std: " << report.psd_std << "\n\n";

    std::cout << "Integrated Charges:\n";
    std::cout << "  Qtot mean: " << std::setprecision(1) << report.qtot_mean << "\n";
    std::cout << "  Qshort mean: " << report.qshort_mean << "\n";
    std::cout << "  Slow component fraction: "
              << (report.slow_component_fraction * 100.0) << "%\n\n";

    std::cout << "Discrimination Performance (Binda Section 6.2):\n";
    std::cout << "  Figure of Merit (FOM): " << std::setprecision(4)
              << report.fom_result.fom << "\n";
    std::cout << "  Reference values:\n";
    std::cout << "    - Binda (NE213 at JET): FOM = 2.25-2.35 (VERY GOOD)\n";
    std::cout << "    - Baselga (CCM baseline): FOM = 0.94\n";
    std::cout << "    - Baselga (best method): FOM = 1.04\n";

    if (report.fom_result.fom > 2.0)
    {
        std::cout << "  [RESULT] VERY GOOD discrimination (matching JET NE213 performance)\n";
        std::cout << "  [RESULT] Estimated neutron events: "
                  << report.estimated_neutron_count << " ("
                  << (report.neutron_fraction * 100.0) << "%)\n";
        std::cout << "  Actual FOM: " << std::setprecision(10)
                  << report.fom_result.fom << "\n";
    }
    else if (report.fom_result.fom > 1.0)
    {
        std::cout << "   [RESULT] GOOD discrimination capability\n";
        std::cout << "  [RESULT] Estimated neutron events: "
                  << report.estimated_neutron_count << " ("
                  << (report.neutron_fraction * 100.0) << "%)\n";
        std::cout << "  Actual FOM: " << std::setprecision(10)
                  << report.fom_result.fom << "\n";
    }
    else if (report.fom_result.fom > 0.7)
    {
        std::cout << "  [RESULT] MODERATE discrimination -> Talk to Felix\n";
        std::cout << "  [RESULT] Estimated neutron events: "
                  << report.estimated_neutron_count << " ("
                  << (report.neutron_fraction * 100.0) << "%)\n";
        std::cout << "  Actual FOM: " << std::setprecision(10)
                  << report.fom_result.fom << "\n";
    }
    else
    {
        std::cout << "  [RESULT] POOR discrimination -> Talk to Felix\n";
        std::cout << "  Actual FOM: " << std::setprecision(10)
                  << report.fom_result.fom << "\n";
    }

    if (report.fom_result.gamma_peak > 0 && report.fom_result.neutron_peak > 0)
    {
        std::cout << "\nPeak Positions:\n";
        std::cout << "  Gamma peak (Y): " << std::setprecision(4)
                  << report.fom_result.gamma_peak << "\n";
        std::cout << "  Neutron peak (n): " << report.fom_result.neutron_peak << "\n";
        std::cout << "  Separation: "
                  << (report.fom_result.neutron_peak - report.fom_result.gamma_peak) << "\n";
        std::cout << "  Peak ratio: " << std::setprecision(2)
                  << (report.fom_result.neutron_peak / report.fom_result.gamma_peak) << "\n";
    }
    std::cout << "\n";
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
