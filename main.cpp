#include "include/WaveformData.h"
#include "include/WaveformProcessor.h"
#include "include/PSDAnalyzer.h"
#include "include/Statistics.h"
#include "include/QtPlotter.h"
#include <iostream>
#include <chrono>
#include <vector>
#include <QApplication>

using namespace ne213;

int main(int argc, char* argv[])
{
    QApplication app(argc, argv);

    std::string filename = "5_mit_waveform_data.txt";
    if (argc > 1)
    {
        filename = argv[1];
    }

    const auto start_time = std::chrono::high_resolution_clock::now();

    try
    {
        std::cout << "Loading waveforms from: " << filename << "\n";
        const auto raw_waveforms = WaveformLoader::load_from_file(filename);
        std::cout << "Loaded " << raw_waveforms.size() << " waveforms from file\n\n";

        std::cout << "Processing waveforms...\n";
        std::vector<WaveformData> processed_waveforms;
        processed_waveforms.reserve(raw_waveforms.size());

        for (const auto& raw : raw_waveforms)
        {
            processed_waveforms.push_back(WaveformProcessor::process_waveform(raw));
        }

        std::cout << "Calculate PSD Parameters...\n";
        const PSDAnalyzer analyzer(25.0, 150.0, 10.0);
        std::vector<PSDParameters> psd_params;
        psd_params.reserve(processed_waveforms.size());

        for (const auto& wf : processed_waveforms)
        {
            psd_params.push_back(analyzer.calculate_parameters(wf));
        }

        std::cout << "Generate statistics report...\n";
        const auto report = Statistics::generate_report(
            processed_waveforms,
            psd_params,
            analyzer.get_short_gate(),
            analyzer.get_total_gate(),
            10.0);

        Statistics::print_report(report);

        const auto fom_result = PSDAnalyzer::calculate_fom(psd_params, 10.0);

        const auto end_time = std::chrono::high_resolution_clock::now();
        const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

        std::cout << "Processing completed in " << duration.count() << "ms\n";

        QtPlotter plotter;

        plotter.plot_raw_waveforms(processed_waveforms);
        plotter.plot_average_pulse_shapes(processed_waveforms, psd_params);
        plotter.plot_psd_scatter(psd_params);
        plotter.plot_psd_histogram(psd_params, 10.0, fom_result);
        plotter.plot_statistics(processed_waveforms, psd_params);

        plotter.show();

        return app.exec();
    }
    catch (std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}