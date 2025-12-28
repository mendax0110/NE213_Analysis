#include "../include/QtPlotter.h"
#include "../include/Statistics.h"
#include "../include/PSDAnalyzer.h"
#include <algorithm>
#include <cmath>
#include <QDateTime>
#include <QValueAxis>
#include <QLineSeries>
#include <QLogValueAxis>
#include <QScatterSeries>
#include <QBarSeries>
#include <QBarCategoryAxis>
#include <QMessageBox>
#include <QStatusBar>
#include <QVBoxLayout>
#include <QMenuBar>
#include <QGroupBox>

using namespace ne213;

QtPlotter::QtPlotter(QWidget *parent) : QMainWindow(parent)
{
    setWindowTitle("NE213 Analysis - Qt Implementation");
    setMinimumSize(1400, 900);

    auto *centralWidget = new QWidget(this);
    setCentralWidget(centralWidget);

    auto mainLayout = new QVBoxLayout(centralWidget);

    tabWidget = new QTabWidget(this);
    mainLayout->addWidget(tabWidget);

    // Tab 1: Raw Waveforms
    auto *rawTab = new QWidget();
    auto *rawLayout = new QHBoxLayout(rawTab);

    rawWaveformChartView = new QChartView();
    rawWaveformChartView->setRenderHint(QPainter::Antialiasing);
    rawLayout->addWidget(rawWaveformChartView, 1);

    avgWaveformChartView = new QChartView();
    avgWaveformChartView->setRenderHint(QPainter::Antialiasing);
    rawLayout->addWidget(avgWaveformChartView, 1);

    tabWidget->addTab(rawTab, "Raw Waveforms");

    // Tab 2: PSD Scatter
    auto *scatterTab = new QWidget();
    auto *scatterLayout = new QVBoxLayout(scatterTab);
    psdScatterChartView = new QChartView();
    psdScatterChartView->setRenderHint(QPainter::Antialiasing);
    scatterLayout->addWidget(psdScatterChartView);
    tabWidget->addTab(scatterTab, "PSD Scatter Plot");

    // Tab 3: Interactive Histogram
    auto *histTab = new QWidget();
    auto *histLayout = new QVBoxLayout(histTab);

    auto *controlGroup = new QGroupBox("PSD Parameters");
    auto *controlLayout = new QGridLayout(controlGroup);

    shortGateSlider = new QSlider(Qt::Horizontal);
    shortGateSlider->setRange(10, 100);
    shortGateSlider->setValue(25);
    shortGateSlider->setTickInterval(10);
    shortGateSlider->setTickPosition(QSlider::TicksBelow);
    controlLayout->addWidget(new QLabel("Short Gate (ns):"), 0, 0);
    controlLayout->addWidget(shortGateSlider, 0, 1);
    shortGateLabel = new QLabel("25.0 ns");
    controlLayout->addWidget(shortGateLabel, 0, 2);

    totalGateSlider = new QSlider(Qt::Horizontal);
    totalGateSlider->setRange(50, 300);
    totalGateSlider->setValue(150);
    totalGateSlider->setTickInterval(25);
    totalGateSlider->setTickPosition(QSlider::TicksBelow);
    controlLayout->addWidget(new QLabel("Total Gate (ns):"), 1, 0);
    controlLayout->addWidget(totalGateSlider, 1, 1);
    totalGateLabel = new QLabel("150.0 ns");
    controlLayout->addWidget(totalGateLabel, 1, 2);

    energyThresholdSlider = new QSlider(Qt::Horizontal);
    energyThresholdSlider->setRange(0, 100);
    energyThresholdSlider->setValue(10);
    energyThresholdSlider->setTickInterval(10);
    energyThresholdSlider->setTickPosition(QSlider::TicksBelow);
    controlLayout->addWidget(new QLabel("Energy Threshold:"), 2, 0);
    controlLayout->addWidget(energyThresholdSlider, 2, 1);
    energyLabel = new QLabel("10.0");
    controlLayout->addWidget(energyLabel, 2, 2);

    binsSlider = new QSlider(Qt::Horizontal);
    binsSlider->setRange(20, 200);
    binsSlider->setValue(100);
    binsSlider->setTickInterval(20);
    binsSlider->setTickPosition(QSlider::TicksBelow);
    controlLayout->addWidget(new QLabel("Histogram Bins:"), 3, 0);
    controlLayout->addWidget(binsSlider, 3, 1);
    binsLabel = new QLabel("100");
    controlLayout->addWidget(binsLabel, 3, 2);

    fomLabel = new QLabel("FOM: --");
    fomLabel->setStyleSheet("font-weight: bold; color: blue;");
    controlLayout->addWidget(fomLabel, 4, 0, 1, 3);

    histLayout->addWidget(controlGroup);

    psdHistogramChartView = new QChartView();
    psdHistogramChartView->setRenderHint(QPainter::Antialiasing);
    histLayout->addWidget(psdHistogramChartView);

    tabWidget->addTab(histTab, "Interactive PSD Histogram");

    // Tab 4: Statistics
    auto *statsTab = new QWidget();
    auto *statsLayout = new QGridLayout(statsTab);

    amplitudeHistView = new QChartView();
    amplitudeHistView->setRenderHint(QPainter::Antialiasing);
    statsLayout->addWidget(amplitudeHistView, 0, 0);

    riseTimeHistView = new QChartView();
    riseTimeHistView->setRenderHint(QPainter::Antialiasing);
    statsLayout->addWidget(riseTimeHistView, 0, 1);

    qtotQshortView = new QChartView();
    qtotQshortView->setRenderHint(QPainter::Antialiasing);
    statsLayout->addWidget(qtotQshortView, 1, 0);

    psdDistView = new QChartView();
    psdDistView->setRenderHint(QPainter::Antialiasing);
    statsLayout->addWidget(psdDistView, 1, 1);

    tabWidget->addTab(statsTab, "Statistical Analysis");

    // Tab 5: Summary
    auto *summaryTab = new QWidget();
    auto *summaryLayout = new QVBoxLayout(summaryTab);
    summaryText = new QTextEdit();
    summaryText->setReadOnly(true);
    summaryText->setFont(QFont("Monospace", 10));
    summaryLayout->addWidget(summaryText);
    tabWidget->addTab(summaryTab, "Summary");

    connect(shortGateSlider, &QSlider::valueChanged, this, &QtPlotter::on_short_gate_changed);
    connect(totalGateSlider, &QSlider::valueChanged, this, &QtPlotter::on_total_gate_changed);
    connect(energyThresholdSlider, &QSlider::valueChanged, this, &QtPlotter::on_energy_threshold_changed);
    connect(binsSlider, &QSlider::valueChanged, this, &QtPlotter::on_bins_changed);

    auto menuBar = new QMenuBar(this);
    setMenuBar(menuBar);

    QMenu *fileMenu = menuBar->addMenu("&File");
    QAction *saveAction = fileMenu->addAction("&Save Plot");
    QAction *exportAction = fileMenu->addAction("&Export Data");
    QAction *exitAction = fileMenu->addAction("&Exit");

    connect(saveAction, &QAction::triggered, this, &QtPlotter::save_plot);
    connect(exportAction, &QAction::triggered, this, &QtPlotter::export_data);
    connect(exitAction, &QAction::triggered, this, &QApplication::quit);

    statusBar()->showMessage("Ready");
}

QtPlotter::~QtPlotter() = default;

std::vector<double> QtPlotter::create_time_axis(const size_t n_samples)
{
    std::vector<double> time(n_samples);
    for (size_t i = 0; i < n_samples; ++i)
    {
        time[i] = i * NS_PRE_SAMPLE;
    }
    return time;
}

QColor QtPlotter::get_color_from_value(const double value, const double min_val, const double max_val)
{
    double norm = (value - min_val) / (max_val - min_val);
    norm = std::max(0.0, std::min(1.0, norm));

    // Blue (low) to Red (high)
    if (norm < 0.5)
    {
        const double t = norm * 2.0;
        return {static_cast<int>(t * 255), static_cast<int>(t * 255), 255};
    }
    else
    {
        const double t = (norm - 0.5) * 2.0;
        return {255, static_cast<int>((1.0 - t) * 255), static_cast<int>((1.0 - t) * 255)};
    }
}

void QtPlotter::plot_average_pulse_shapes(const std::vector<WaveformData>& waveforms, const std::vector<PSDParameters>& psd_params) const
{
    if (waveforms.empty() || psd_params.empty()) return;

    std::vector<std::pair<double, size_t>> psd_indices;
    for (size_t i = 0; i < psd_params.size(); ++i)
    {
        if (psd_params[i].is_valid)
        {
            psd_indices.emplace_back(psd_params[i].psd_values, i);
        }
    }

    std::ranges::sort(psd_indices);
    if (psd_indices.size() < 10) return;
    const size_t third = psd_indices.size() / 3;

    std::vector avg_low(waveforms[0].data.size(), 0.0);
    std::vector avg_high(waveforms[0].data.size(), 0.0);

    for (size_t i = 0; i < third; ++i)
    {
        const size_t idx = psd_indices[i].second;
        for (size_t j = 0; j < avg_low.size(); ++j)
        {
            avg_low[j] += waveforms[idx].data[j];
        }
    }

    for (size_t i = psd_indices.size() - third; i < psd_indices.size(); ++i)
    {
        const size_t idx = psd_indices[i].second;
        for (size_t j = 0; j < avg_high.size(); ++j)
        {
            avg_high[j] += waveforms[idx].data[j];
        }
    }

    for (auto& v : avg_low) v /= third;
    for (auto& v : avg_high) v /= third;

    const auto time_axis = create_time_axis(avg_low.size());

    const auto chart = new QChart();
    chart->setTitle("Average Pulse Shapes: Neutron vs Gamma (log scale)");
    chart->setAnimationOptions(QChart::NoAnimation);
    chart->setMargins(QMargins(10, 10, 10, 10));

    const auto lowSeries = new QLineSeries();
    lowSeries->setName("Low PSD (Gamma)");
    lowSeries->setPen(QPen(Qt::blue, 2.0));

    for (size_t i = 0; i < avg_low.size(); ++i)
    {
        lowSeries->append(time_axis[i], std::max(0.1, avg_low[i]));
    }

    const auto highSeries = new QLineSeries();
    highSeries->setName("High PSD (Neutron)");
    highSeries->setPen(QPen(Qt::red, 2.0));

    for (size_t i = 0; i < avg_high.size(); ++i)
    {
        highSeries->append(time_axis[i], std::max(0.1, avg_high[i]));
    }

    chart->addSeries(lowSeries);
    chart->addSeries(highSeries);

    add_gate_lines(chart, currentShortGate, currentTotalGate, 100000);

    const auto axisX = new QValueAxis();
    axisX->setTitleText("Time (ns)");
    axisX->setLabelFormat("%.0f");
    axisX->setRange(0, 175);
    axisX->setTickCount(8);
    chart->addAxis(axisX, Qt::AlignBottom);

    const auto axisY = new QLogValueAxis();
    axisY->setTitleText("Amplitude");
    axisY->setLabelFormat("%.0e");
    axisY->setBase(10.0);
    axisY->setMin(0.1);
    axisY->setMax(100000);
    chart->addAxis(axisY, Qt::AlignLeft);

    lowSeries->attachAxis(axisX);
    lowSeries->attachAxis(axisY);
    highSeries->attachAxis(axisX);
    highSeries->attachAxis(axisY);

    chart->legend()->setVisible(true);
    chart->legend()->setAlignment(Qt::AlignBottom);
    avgWaveformChartView->setChart(chart);
}

void QtPlotter::plot_raw_waveforms(const std::vector<WaveformData>& waveforms, const double short_gate_ns, const double total_gate_ns)
{
    if (waveforms.empty()) return;

    currentWaveforms = waveforms;
    currentShortGate = short_gate_ns;
    currentTotalGate = total_gate_ns;

    const auto rawChart = new QChart();
    rawChart->setTitle("NE213 Scintillation Pulses (D-D Fusion)");
    rawChart->setAnimationOptions(QChart::NoAnimation);
    rawChart->setMargins(QMargins(10, 10, 10, 10));

    const auto time_axis = create_time_axis(waveforms[0].data.size());

    for (const auto & wf : waveforms)
    {
        const auto series = new QLineSeries();
        series->setPen(QPen(QColor(0, 0, 255, 77), 0.5));

        for (size_t j = 0; j < wf.data.size(); ++j)
        {
            series->append(time_axis[j], wf.data[j]);
        }
        rawChart->addSeries(series);
    }

    add_gate_lines(rawChart, short_gate_ns, total_gate_ns, 25000);

    const auto axisX = new QValueAxis();
    axisX->setTitleText("Time (ns)");
    axisX->setLabelFormat("%.0f");
    axisX->setRange(0, 175);
    axisX->setTickCount(8);
    rawChart->addAxis(axisX, Qt::AlignBottom);

    const auto axisY = new QValueAxis();
    axisY->setTitleText("Amplitude (baseline corrected)");
    axisY->setLabelFormat("%.0f");
    axisY->setRange(0, 25000);
    axisY->setTickCount(6);
    rawChart->addAxis(axisY, Qt::AlignLeft);

    for (const auto series : rawChart->series())
    {
        series->attachAxis(axisX);
        series->attachAxis(axisY);
    }

    rawChart->legend()->setVisible(true);
    rawChart->legend()->setAlignment(Qt::AlignBottom);

    rawWaveformChartView->setChart(rawChart);
}

void QtPlotter::plot_psd_scatter(const std::vector<PSDParameters>& psd_params)
{
    currentPsdParams = psd_params;

    std::vector<double> qtot, psd_values;
    for (const auto& p : psd_params)
    {
        if (p.is_valid)
        {
            qtot.push_back(p.qtot);
            psd_values.push_back(p.psd_values);
        }
    }

    if (qtot.empty()) return;

    const auto chart = new QChart();
    chart->setTitle("NE213 Pulse Shape Discrimination (Binda Method)\nD-D Fusion\n(Neutrons: high PSD, Gammas: low PSD)");
    chart->setAnimationOptions(QChart::NoAnimation);

    const double min_psd = *std::ranges::min_element(psd_values);
    const double max_psd = *std::ranges::max_element(psd_values);

    constexpr int num_color_bands = 20;
    std::vector<QScatterSeries*> series_list(num_color_bands);

    for (int s = 0; s < num_color_bands; ++s)
    {
        series_list[s] = new QScatterSeries();
        const double norm = static_cast<double>(s) / (num_color_bands - 1);
        QColor color = get_color_from_value(norm * (max_psd - min_psd) + min_psd, min_psd, max_psd);
        series_list[s]->setColor(color);
        series_list[s]->setBorderColor(Qt::black);
        series_list[s]->setMarkerSize(6.0);

        const double psd_low = min_psd + (s * (max_psd - min_psd) / num_color_bands);
        const double psd_high = min_psd + ((s + 1) * (max_psd - min_psd) / num_color_bands);
        series_list[s]->setName(QString("%1-%2").arg(psd_low, 0, 'f', 3).arg(psd_high, 0, 'f', 3));
    }

    for (size_t i = 0; i < qtot.size(); ++i)
    {
        const double norm = (psd_values[i] - min_psd) / (max_psd - min_psd);
        const int band_idx = std::min(num_color_bands - 1, static_cast<int>(norm * num_color_bands));
        series_list[band_idx]->append(qtot[i], psd_values[i]);
    }

    for (auto* s : series_list)
    {
        if (s->count() > 0)
        {
            chart->addSeries(s);
        }
    }

    std::vector<double> sorted_psd = psd_values;
    std::ranges::sort(sorted_psd);
    const double median_psd = sorted_psd[sorted_psd.size() / 2];

    const double min_qtot = *std::ranges::min_element(qtot);
    const double max_qtot = *std::ranges::max_element(qtot);

    const auto medianLine = new QLineSeries();
    medianLine->setName(QString("Discrimination threshold (median PSD=%1)").arg(median_psd, 0, 'f', 3));
    medianLine->append(min_qtot, median_psd);
    medianLine->append(max_qtot, median_psd);
    medianLine->setPen(QPen(Qt::green, 2.0, Qt::DashLine));

    chart->addSeries(medianLine);

    const auto axisX = new QValueAxis();
    axisX->setTitleText("Qtot - Total Charge (Energy Proxy)");
    axisX->setLabelFormat("%.0f");
    axisX->setRange(100000, 700000);
    axisX->setTickCount(7);
    chart->addAxis(axisX, Qt::AlignBottom);

    const auto axisY = new QValueAxis();
    axisY->setTitleText("PSD = (Qtot - Qshort) / Qtot");
    axisY->setLabelFormat("%.3f");
    axisY->setRange(0.125, 0.300);
    axisY->setTickCount(8);
    chart->addAxis(axisY, Qt::AlignLeft);

    for (auto* s : series_list)
    {
        if (s->count() > 0)
        {
            s->attachAxis(axisX);
            s->attachAxis(axisY);
        }
    }
    medianLine->attachAxis(axisX);
    medianLine->attachAxis(axisY);

    chart->legend()->setVisible(true);
    chart->legend()->setAlignment(Qt::AlignRight);

    psdScatterChartView->setChart(chart);
}

void QtPlotter::plot_psd_histogram(const std::vector<PSDParameters>& psd_params,
                                 const double energy_threshold,
                                 const FOMResult& fom_result,
                                 const double short_gate_ns,
                                 const double total_gate_ns)
{
    currentPsdParams = psd_params;
    currentEnergyThreshold = energy_threshold;
    currentFomResult = fom_result;
    currentShortGate = short_gate_ns;
    currentTotalGate = total_gate_ns;

    update_psd_histogram();
}

void QtPlotter::update_psd_histogram()
{
    const PSDAnalyzer analyzer(currentShortGate, currentTotalGate, 10.0);

    std::vector<PSDParameters> updated_params;
    for (const auto& wf : currentWaveforms)
    {
        updated_params.push_back(analyzer.calculate_parameters(wf));
    }

    currentFomResult = PSDAnalyzer::calculate_fom(updated_params, currentEnergyThreshold);

    std::vector<double> filtered_psd;
    for (const auto& p : updated_params)
    {
        if (p.is_valid && p.qtot > currentEnergyThreshold)
        {
            filtered_psd.push_back(p.psd_values);
        }
    }

    if (filtered_psd.empty()) return;

    const double min_psd = *std::ranges::min_element(filtered_psd);
    const double max_psd = *std::ranges::max_element(filtered_psd);
    const double bin_width = (max_psd - min_psd) / currentBins;

    std::vector counts(currentBins, 0);
    std::vector<double> bin_centers(currentBins);

    for (int i = 0; i < currentBins; ++i)
    {
        bin_centers[i] = min_psd + (i + 0.5) * bin_width;
    }

    for (const double psd_val : filtered_psd)
    {
        int bin_idx = static_cast<int>((psd_val - min_psd) / bin_width);
        if (bin_idx >= currentBins) bin_idx = currentBins - 1;
        if (bin_idx >= 0 && bin_idx < currentBins)
        {
            counts[bin_idx]++;
        }
    }

    auto *chart = new QChart();

    const int short_gate_samples = ns_to_sample(currentShortGate);
    const int total_gate_samples = ns_to_sample(currentTotalGate);

    QString title = QString("NE213 PSD Distribution (FOM=%1)\n"
                           "Valid events: %2 | Short gate: %3 ns (%4 samples) | Total gate: %5 ns (%6 samples)")
                    .arg(currentFomResult.fom, 0, 'f', 3)
                    .arg(filtered_psd.size())
                    .arg(currentShortGate, 0, 'f', 1)
                    .arg(short_gate_samples)
                    .arg(currentTotalGate, 0, 'f', 1)
                    .arg(total_gate_samples);

    if (currentFomResult.fom > 0)
    {
        title += QString("\nGamma peak=%1 | Neutron peak=%2")
                 .arg(currentFomResult.gamma_peak, 0, 'f', 3)
                 .arg(currentFomResult.neutron_peak, 0, 'f', 3);

        if (currentFomResult.fom > 2.0)
        {
            title += " | VERY GOOD (matching Binda JET results)";
        }
        else if (currentFomResult.fom > 1.0)
        {
            title += " | GOOD discrimination";
        }
        else if (currentFomResult.fom > 0.5)
        {
            title += " | MODERATE discrimination";
        }
        else
        {
            title += " | POOR - optimize windows";
        }
    }

    chart->setTitle(title);
    chart->setAnimationOptions(QChart::NoAnimation);

    auto series = new QBarSeries();
    auto barSet = new QBarSet("PSD Distribution");

    for (int i = 0; i < currentBins; ++i)
    {
        *barSet << counts[i];
    }

    barSet->setColor(QColor(0, 0, 255, 180));
    barSet->setBorderColor(Qt::black);
    series->append(barSet);

    chart->addSeries(series);

    QStringList categories;
    const int label_step = std::max(1, currentBins / 10);

    for (int i = 0; i < currentBins; ++i)
    {
        if (i % label_step == 0)
        {
            categories << QString::number(bin_centers[i], 'f', 3);
        }
        else
        {
            categories << "";
        }
    }

    auto categoryAxis = new QBarCategoryAxis();
    categoryAxis->append(categories);
    categoryAxis->setTitleText("PSD = (Qtot - Qshort) / Qtot");
    chart->addAxis(categoryAxis, Qt::AlignBottom);
    series->attachAxis(categoryAxis);

    auto axisY = new QValueAxis();
    axisY->setTitleText("Counts");
    axisY->setLabelFormat("%d");
    chart->addAxis(axisY, Qt::AlignLeft);
    series->attachAxis(axisY);

    chart->legend()->setVisible(false);

    fomLabel->setText(QString("Figure of Merit (FOM): %1").arg(currentFomResult.fom, 0, 'f', 4));

    psdHistogramChartView->setChart(chart);
}

void QtPlotter::plot_statistics(const std::vector<WaveformData>& waveforms, const std::vector<PSDParameters>& psd_params) const
{
    if (waveforms.empty() || psd_params.empty()) return;

    // 1. Amplitude distribution
    std::vector<double> amplitudes;
    for (const auto& wf : waveforms)
    {
        amplitudes.push_back(wf.max_amplitude);
    }

    QChart *ampChart = create_value_histogram(amplitudes, 50,
                                            "Pulse Amplitude Distribution",
                                            "Peak Amplitude (a.u.)",
                                            "Counts",
                                            QColor(128, 0, 128));

    auto ampXAxes = ampChart->axes(Qt::Horizontal);
    if (!ampXAxes.isEmpty())
    {
        auto axisX = qobject_cast<QBarCategoryAxis*>(ampXAxes.first());
        if (axisX)
        {
            QStringList ampLabels;
            for (int val = 0; val <= 25000; val += 5000)
            {
                ampLabels << QString::number(val);
            }
            axisX->clear();
            axisX->append(ampLabels);
        }
    }

    auto ampYAxes = ampChart->axes(Qt::Vertical);
    if (!ampYAxes.isEmpty())
    {
        auto axisY = qobject_cast<QValueAxis*>(ampYAxes.first());
        if (axisY)
        {
            axisY->setRange(0, 800);
            axisY->setTickCount(9);
        }
    }

    amplitudeHistView->setChart(ampChart);

    // 2. Rise time distribution
    std::vector<double> rise_times_ns;
    for (const auto& wf : waveforms)
    {
        rise_times_ns.push_back(samples_to_ns(wf.rise_time));
    }

    QChart *riseChart = create_value_histogram(rise_times_ns, 50,
                                             "Pulse Rise Time Distribution\n(10%-90% of peak)",
                                             "Rise Time (ns)",
                                             "Counts",
                                             QColor(255, 165, 0));

    auto riseXAxes = riseChart->axes(Qt::Horizontal);
    if (!riseXAxes.isEmpty())
    {
        auto axisX = qobject_cast<QBarCategoryAxis*>(riseXAxes.first());
        if (axisX)
        {
            QStringList riseLabels;
            for (int val = 0; val <= 175; val += 25)
            {
                riseLabels << QString::number(val);
            }
            axisX->clear();
            axisX->append(riseLabels);
        }
    }

    auto riseYAxes = riseChart->axes(Qt::Vertical);
    if (!riseYAxes.isEmpty())
    {
        if (auto axisY = qobject_cast<QValueAxis*>(riseYAxes.first()))
        {
            axisY->setRange(0, 1400);
            axisY->setTickCount(8);
        }
    }

    riseTimeHistView->setChart(riseChart);

    // 3. Qtot vs Qshort scatter
    std::vector<double> qtot_plot, qshort_plot, psd_plot;
    for (const auto&[psd_values, qtot, qshort, is_valid] : psd_params)
    {
        if (is_valid)
        {
            qtot_plot.push_back(qtot);
            qshort_plot.push_back(qshort);
            psd_plot.push_back(psd_values);
        }
    }

    if (!qtot_plot.empty())
    {
        auto scatterChart = new QChart();
        scatterChart->setTitle("Short vs Total Charge (Binda Method)\n(Valid events, pile-up rejected)");
        scatterChart->setAnimationOptions(QChart::NoAnimation);

        const double min_psd = *std::ranges::min_element(psd_plot);
        const double max_psd = *std::ranges::max_element(psd_plot);

        constexpr int num_color_bands = 20;
        std::vector<QScatterSeries*> series_list(num_color_bands);

        for (int s = 0; s < num_color_bands; ++s)
        {
            series_list[s] = new QScatterSeries();
            double norm = static_cast<double>(s) / (num_color_bands - 1);
            QColor color = get_color_from_value(norm * (max_psd - min_psd) + min_psd, min_psd, max_psd);
            series_list[s]->setColor(color);
            series_list[s]->setBorderColor(Qt::black);
            series_list[s]->setMarkerSize(4.0);
        }

        for (size_t i = 0; i < qtot_plot.size(); ++i)
        {
            double norm = (psd_plot[i] - min_psd) / (max_psd - min_psd);
            int band_idx = std::min(num_color_bands - 1, static_cast<int>(norm * num_color_bands));
            series_list[band_idx]->append(qtot_plot[i], qshort_plot[i]);
        }

        for (auto* s : series_list)
        {
            if (s->count() > 0)
            {
                scatterChart->addSeries(s);
            }
        }

        auto axisX = new QValueAxis();
        axisX->setTitleText("Qtot (total charge)");
        axisX->setLabelFormat("%.0f");
        axisX->setRange(100000, 700000);
        axisX->setTickCount(7);  // 100k, 200k, ..., 700k
        scatterChart->addAxis(axisX, Qt::AlignBottom);

        auto axisY = new QValueAxis();
        axisY->setTitleText("Qshort (short gate charge)");
        axisY->setLabelFormat("%.0f");
        axisY->setRange(100000, 500000);
        axisY->setTickCount(5);  // 100k, 200k, 300k, 400k, 500k
        scatterChart->addAxis(axisY, Qt::AlignLeft);

        for (auto* s : series_list)
        {
            if (s->count() > 0)
            {
                s->attachAxis(axisX);
                s->attachAxis(axisY);
            }
        }

        scatterChart->legend()->setVisible(false);
        qtotQshortView->setChart(scatterChart);
    }

    // 4. PSD parameter distribution (log scale)
    std::vector<double> psd_valid;
    for (const auto& p : psd_params)
    {
        if (p.is_valid)
        {
            psd_valid.push_back(p.psd_values);
        }
    }

    if (!psd_valid.empty())
    {
        constexpr int bins = 100;
        const double min_psd = *std::ranges::min_element(psd_valid);
        const double max_psd = *std::ranges::max_element(psd_valid);

        std::vector<double> bin_edges(bins + 1);
        const double range = max_psd - min_psd;

        for (int i = 0; i <= bins; ++i)
        {
            bin_edges[i] = min_psd + (i * range / bins);
        }

        std::vector<int> counts(bins, 0);
        std::vector<double> bin_centers(bins);


        for (int i = 0; i < bins; ++i)
        {
            bin_centers[i] = (bin_edges[i] + bin_edges[i + 1]) / 2.0;
        }

        for (const double psd_val : psd_valid)
        {
            if (psd_val < min_psd || psd_val > max_psd) continue;

            int bin_idx = static_cast<int>((psd_val - min_psd) / range * bins);

            if (bin_idx >= bins)
            {
                bin_idx = bins - 1;
            }

            counts[bin_idx]++;
        }

        auto psdChart = new QChart();
        psdChart->setTitle("PSD Parameter Distribution (log scale)");
        psdChart->setAnimationOptions(QChart::NoAnimation);

        auto series = new QBarSeries();
        auto barSet = new QBarSet("PSD");

        for (int i = 0; i < bins; ++i)
        {
            *barSet << std::max(1, counts[i]);
        }

        barSet->setColor(QColor(0, 128, 0, 180));
        barSet->setBorderColor(Qt::black);
        series->append(barSet);
        psdChart->addSeries(series);


        QStringList psdLabels;
        constexpr int label_step = bins / 8;
        for (int i = 0; i < bins; ++i)
        {
            if (i % label_step == 0)
            {
                psdLabels << QString::number(bin_centers[i], 'f', 3);
            }
            else
            {
                psdLabels << "";
            }
        }

        auto *axisX = new QBarCategoryAxis();
        axisX->append(psdLabels);
        axisX->setTitleText("PSD = (Qtot - Qshort) / Qtot");
        psdChart->addAxis(axisX, Qt::AlignBottom);
        series->attachAxis(axisX);

        // Log Y axis
        auto *logAxis = new QLogValueAxis();
        logAxis->setTitleText("Counts");
        logAxis->setBase(10.0);
        logAxis->setLabelFormat("%.0e");
        logAxis->setMin(1);      // 10⁰
        logAxis->setMax(10);     // 10¹
        psdChart->addAxis(logAxis, Qt::AlignLeft);
        series->attachAxis(logAxis);

        psdChart->legend()->setVisible(false);

        psdDistView->setChart(psdChart);
    }
}

QChart* QtPlotter::create_value_histogram(const std::vector<double>& data,
                                        const int bins,
                                        const QString& title,
                                        const QString& xlabel,
                                        const QString& ylabel,
                                        const QColor& color)
{
    if (data.empty()) return nullptr;

    const double min_val = *std::ranges::min_element(data);
    const double max_val = *std::ranges::max_element(data);
    const double bin_width = (max_val - min_val) / bins;

    std::vector counts(bins, 0);
    std::vector<double> bin_centers(bins);

    for (int i = 0; i < bins; ++i)
    {
        bin_centers[i] = min_val + (i + 0.5) * bin_width;
    }

    for (const double val : data)
    {
        const int bin_idx = static_cast<int>((val - min_val) / bin_width);
        if (bin_idx >= 0 && bin_idx < bins)
        {
            counts[bin_idx]++;
        }
    }

    const auto chart = new QChart();
    chart->setTitle(title);
    chart->setAnimationOptions(QChart::NoAnimation);

    // Create BAR series like Python matplotlib
    const auto series = new QBarSeries();
    const auto barSet = new QBarSet("Data");

    for (int i = 0; i < bins; ++i)
    {
        *barSet << counts[i];
    }

    barSet->setColor(color);
    barSet->setBorderColor(Qt::black);
    series->append(barSet);
    chart->addSeries(series);

    // Create X axis categories - show every Nth label to avoid clutter
    QStringList categories;
    const int label_step = std::max(1, bins / 8);  // Show ~8 labels

    for (int i = 0; i < bins; ++i)
    {
        if (i % label_step == 0)
        {
            categories << QString::number(bin_centers[i], 'f', 0);
        }
        else
        {
            categories << "";  // Empty label for bars in between
        }
    }

    auto *axisX = new QBarCategoryAxis();
    axisX->append(categories);
    axisX->setTitleText(xlabel);
    chart->addAxis(axisX, Qt::AlignBottom);
    series->attachAxis(axisX);

    // Y axis
    auto *axisY = new QValueAxis();
    axisY->setTitleText(ylabel);
    axisY->setLabelFormat("%d");
    chart->addAxis(axisY, Qt::AlignLeft);
    series->attachAxis(axisY);

    chart->legend()->setVisible(false);

    return chart;
}

void QtPlotter::add_gate_lines(QChart* chart, const double short_gate_ns, const double total_gate_ns, const double max_y)
{
    const auto shortSeries = new QLineSeries();
    shortSeries->setName(QString("Short gate end (%1 ns)").arg(short_gate_ns));
    shortSeries->append(short_gate_ns, 0);
    shortSeries->append(short_gate_ns, max_y);
    shortSeries->setPen(QPen(QColor(0, 128, 0, 128), 1.5, Qt::DashLine));
    chart->addSeries(shortSeries);

    const auto totalSeries = new QLineSeries();
    totalSeries->setName(QString("Total gate end (%1 ns)").arg(total_gate_ns));
    totalSeries->append(total_gate_ns, 0);
    totalSeries->append(total_gate_ns, max_y);
    totalSeries->setPen(QPen(QColor(128, 0, 0, 128), 1.5, Qt::DashLine));
    chart->addSeries(totalSeries);
}

void QtPlotter::on_short_gate_changed(const int value)
{
    currentShortGate = static_cast<double>(value);
    shortGateLabel->setText(QString("%1 ns").arg(value));
    update_psd_histogram();
}

void QtPlotter::on_total_gate_changed(const int value)
{
    currentTotalGate = static_cast<double>(value);
    totalGateLabel->setText(QString("%1 ns").arg(value));
    update_psd_histogram();
}

void QtPlotter::on_energy_threshold_changed(const int value)
{
    currentEnergyThreshold = static_cast<double>(value);
    energyLabel->setText(QString::number(value));
    update_psd_histogram();
}

void QtPlotter::on_bins_changed(const int value)
{
    currentBins = value;
    binsLabel->setText(QString::number(value));
    update_psd_histogram();
}

void QtPlotter::save_plot()
{
    const QString fileName = QFileDialog::getSaveFileName(this,
        "Save Plot", "", "PNG Files (*.png);;PDF Files (*.pdf);;SVG Files (*.svg)");

    if (!fileName.isEmpty())
    {
        QChartView *currentView = nullptr;

        const int currentTab = tabWidget->currentIndex();
        switch (currentTab)
        {
            case 0: currentView = rawWaveformChartView; break;
            case 1: currentView = psdScatterChartView; break;
            case 2: currentView = psdHistogramChartView; break;
            default: break;
        }

        if (currentView && currentView->chart())
        {
            const QPixmap pixmap = currentView->grab();
            if (pixmap.save(fileName))
            {
                statusBar()->showMessage("Plot saved successfully", 3000);
            }
            else
            {
                QMessageBox::warning(this, "Error", "Failed to save plot");
            }
        }
    }
}

void QtPlotter::export_data()
{
    const QString fileName = QFileDialog::getSaveFileName(this,
        "Export Data", "", "CSV Files (*.csv);;Text Files (*.txt)");

    if (!fileName.isEmpty())
    {
        QFile file(fileName);
        if (file.open(QIODevice::WriteOnly | QIODevice::Text))
        {
            QTextStream stream(&file);

            stream << "Event,Amplitude,RiseTime,PSD,Qtot,Qshort,Valid\n";

            for (size_t i = 0; i < currentPsdParams.size(); ++i)
            {
                const auto&[psd_values, qtot, qshort, is_valid] = currentPsdParams[i];
                stream << i << ","
                       << (i < currentWaveforms.size() ? currentWaveforms[i].max_amplitude : 0) << ","
                       << (i < currentWaveforms.size() ? samples_to_ns(currentWaveforms[i].rise_time) : 0) << ","
                       << psd_values << ","
                       << qtot << ","
                       << qshort << ","
                       << (is_valid ? "Yes" : "No") << "\n";
            }

            file.close();
            statusBar()->showMessage("Data exported successfully", 3000);
        }
        else
        {
            QMessageBox::warning(this, "Error", "Failed to export data");
        }
    }
}

void QtPlotter::show_summary(const std::vector<WaveformData>& waveforms, const std::vector<PSDParameters>& psd_params) const
{
    const StatisticsReport report = Statistics::generate_report(waveforms, psd_params, currentShortGate, currentTotalGate, currentEnergyThreshold);

    const std::stringstream ss;
    Statistics::print_report(report);

    summaryText->setText(QString::fromStdString(ss.str()));
}