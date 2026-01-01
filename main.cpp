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

    std::string filename;

    if (argc > 1)
    {
        filename = argv[1];
    }
    else
    {
        QString selected = QtPlotter::select_file();
        if (selected.isEmpty())
        {
            std::cerr << "No file selected. Exiting.\n";
            return 1;
        }
        filename = selected.toStdString();
    }

    const auto start_time = std::chrono::high_resolution_clock::now();

    try
    {
        std::cout << "Loading waveforms from: " << filename << "\n";
        const auto raw_waveforms = WaveformLoader::load_from_file(filename);
        std::cout << "Loaded " << raw_waveforms.size() << " waveforms from file\n";

        const FileType file_type = WaveformLoader::classify_file_type(filename);
        std::cout << "File classification: " << WaveformLoader::file_type_to_string(file_type) << "\n\n";

        std::cout << "Processing waveforms...\n";
        std::vector<WaveformData> processed_waveforms;
        processed_waveforms.reserve(raw_waveforms.size());

        for (const auto& raw : raw_waveforms)
        {
            processed_waveforms.push_back(WaveformProcessor::process_waveform(raw));
        }

        std::vector<PSDParameters> psd_params;
        FOMResult fom_result;

        if (file_type == FileType::MIT)
        {
            std::cout << "Calculating PSD parameters...\n";
            const PSDAnalyzer analyzer(25.0, 150.0, 10.0);
            psd_params.reserve(processed_waveforms.size());

            for (const auto& wf : processed_waveforms)
            {
                psd_params.push_back(analyzer.calculate_parameters(wf));
            }

            fom_result = PSDAnalyzer::calculate_fom(psd_params, 10.0);

            std::cout << "Generating statistics report...\n";
            const auto report = Statistics::generate_report(
                processed_waveforms,
                psd_params,
                analyzer.get_short_gate(),
                analyzer.get_total_gate(),
                10.0,
                file_type
            );

            Statistics::print_report(report);
        }
        else
        {
            std::cout << "Background file detected - skipping PSD calculations\n\n";

            psd_params.resize(processed_waveforms.size());
            for (auto& p : psd_params)
            {
                p.is_valid = true;
            }

            const auto report = Statistics::generate_report(
                processed_waveforms,
                psd_params,
                25.0,
                150.0,
                10.0,
                file_type
            );

            Statistics::print_report(report);
        }

        const auto end_time = std::chrono::high_resolution_clock::now();
        const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

        std::cout << "Processing completed in " << duration.count() << "ms\n";

        QtPlotter plotter;
        plotter.set_file_type(file_type);

        plotter.plot_raw_waveforms(processed_waveforms);
        plotter.plot_average_pulse_shapes(processed_waveforms, psd_params);
        plotter.plot_psd_scatter(psd_params);
        plotter.plot_psd_histogram(psd_params, 10.0, fom_result);
        plotter.plot_statistics(processed_waveforms, psd_params);
        plotter.show_summary(processed_waveforms, psd_params);

        plotter.show();

        return app.exec();
    }
    catch (std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}