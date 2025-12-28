#pragma once

#include <QMainWindow>
#include <QTabWidget>
#include <QWidget>
#include <QChartView>
#include <QBarSet>
#include <QLabel>
#include <QApplication>
#include <QListWidget>
#include <QTextEdit>
#include <QFileDialog>
#include <QAreaSeries>
#include <QBrush>
#include <QColor>
#include <cmath>
#include <vector>
#include <memory>
#include "WaveformData.h"

QT_USE_NAMESPACE

/// @brief Namespace for NE213 detector data analysis \namespace ne213
namespace ne213
{
    /// @brief
    class QtPlotter : public QMainWindow
    {
        Q_OBJECT
    public:
        /**
         * @brief Constructor
         * @param parent The parent widget
         */
        explicit QtPlotter(QWidget *parent = nullptr);

        /**
         * @brief Destructor
         */
        ~QtPlotter();

        /**
         * @brief Plot raw waveforms
         * @param waveforms The processed waveforms
         * @param short_gate_ns The short gate duration in nanoseconds
         * @param total_gate_ns The total gate duration in nanoseconds
         */
        void plot_raw_waveforms(const std::vector<WaveformData>& waveforms, double short_gate_ns = 25.0, double total_gate_ns = 150.0);

        /**
         * @brief Plot average pulse shapes
         * @param waveforms The processed waveforms
         * @param psd_params The corresponding PSD parameters
         */
        void plot_average_pulse_shapes(const std::vector<WaveformData>& waveforms, const std::vector<PSDParameters>& psd_params) const;

        /**
         * @brief Plot PSD scatter plot
         * @param psd_params The PSD parameters to plot
         */
        void plot_psd_scatter(const std::vector<PSDParameters>& psd_params);

        /**
         * @brief Plot interactive PSD histogram with adjustable parameters
         * @param psd_params The PSD parameters to plot
         * @param energy_threshold THe energy threshold for event selection
         * @param fom_result The Figure of Merit result
         * @param short_gate_ns THe short gate duration in nanoseconds
         * @param total_gate_ns The total gate duration in nanoseconds
         */
        void plot_psd_histogram(const std::vector<PSDParameters>& psd_params,
                                double energy_threshold = 10.0,
                                const FOMResult& fom_result = FOMResult(),
                                double short_gate_ns = 25.0,
                                double total_gate_ns = 150.0);

        /**
         * @brief Plot various statistical histograms
         * @param waveforms The processed waveforms
         * @param psd_params The corresponding PSD parameters
         */
        void plot_statistics(const std::vector<WaveformData>& waveforms, const std::vector<PSDParameters>& psd_params) const;

        /**
         * @brief Show summary report
         * @param waveforms The processed waveforms
         * @param psd_params The corresponding PSD parameters
         */
        void show_summary(const std::vector<WaveformData>& waveforms, const std::vector<PSDParameters>& psd_params) const;

    private slots:
        void on_short_gate_changed(int value);
        void on_total_gate_changed(int value);
        void on_energy_threshold_changed(int value);
        void on_bins_changed(int value);
        void update_psd_histogram();
        void save_plot();
        void export_data();

    private:
        // UI Components
        QTabWidget *tabWidget;

        // Raw waveforms tab
        QChartView *rawWaveformChartView;
        QChartView *avgWaveformChartView;

        // PSD scatter tab
        QChartView *psdScatterChartView;

        // Interactive histogram tab
        QChartView *psdHistogramChartView;
        QSlider *shortGateSlider;
        QSlider *totalGateSlider;
        QSlider *energyThresholdSlider;
        QSlider *binsSlider;
        QLabel *shortGateLabel;
        QLabel *totalGateLabel;
        QLabel *energyLabel;
        QLabel *binsLabel;
        QLabel *fomLabel;

        // Statistics tab
        QChartView *amplitudeHistView;
        QChartView *riseTimeHistView;
        QChartView *qtotQshortView;
        QChartView *psdDistView;

        // Summary tab
        QTextEdit *summaryText;

        // Data storage
        std::vector<WaveformData> currentWaveforms;
        std::vector<PSDParameters> currentPsdParams;
        double currentShortGate = 25.0;
        double currentTotalGate = 150.0;
        double currentEnergyThreshold = 10.0;
        int currentBins = 100;
        FOMResult currentFomResult;

        /**
         * @brief Create time axis in nanoseconds
         * @param n_samples Number of samples
         * @return Vector of time values in nanoseconds
         */
        static std::vector<double> create_time_axis(size_t n_samples);

        /**
         * @brief Map a value to a color gradient between blue and red
         * @param value The value to map
         * @param min_val Minimum value of the range
         * @param max_val Maximum value of the range
         * @return Corresponding QColor
         */
        static QColor get_color_from_value(double value, double min_val, double max_val);

        /**
         * @brief Create a histogram chart from data
         * @param data Input data for histogram
         * @param bins Number of bins
         * @param title Chart title
         * @param xlabel X-axis label
         * @param ylabel Y-axis label
         * @param color Bar color
         * @return Pointer to the created QChart
         */
        static QChart* create_value_histogram(const std::vector<double>& data,
                                              int bins,
                                              const QString& title,
                                              const QString& xlabel,
                                              const QString& ylabel,
                                              const QColor& color);

        /**
         * @brief Add gate lines to the chart
         * @param chart The chart to add lines to
         * @param short_gate_ns The short gate duration in nanoseconds
         * @param total_gate_ns The total gate duration in nanoseconds
         * @param max_y The maximum y-value for the lines
         */
        static void add_gate_lines(QChart* chart, double short_gate_ns, double total_gate_ns, double max_y);
    };
}