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
#include <QSlider>
#include <cmath>
#include <vector>
#include <memory>
#include "WaveformData.h"

QT_USE_NAMESPACE

/// @brief Namespace for NE213 detector data analysis \namespace ne213
namespace ne213
{
    /**
     * @brief Main window for NE213 analysis visualization
     *
     * Provides a tabbed interface for viewing:
     * - Raw waveforms and average pulse shapes
     * - 2D PSD scatter plot with colormap
     * - Interactive PSD histogram with parameter sliders
     * - Statistical distribution plots
     * - Analysis summary report
     *
     * Supports file selection, plot saving, and data export.
     */
    class QtPlotter : public QMainWindow
    {
        Q_OBJECT

    public:
        /**
         * @brief Construct the main plotter window
         * @param parent Parent widget (nullptr for top-level window)
         */
        explicit QtPlotter(QWidget* parent = nullptr);

        /**
         * @brief Destructor
         */
        ~QtPlotter() override;

        /**
         * @brief Open file selection dialog
         *
         * Opens a QFileDialog for selecting waveform data files.
         * Returns the selected filename or empty string if cancelled.
         *
         * @return Selected filename, or empty string if cancelled
         */
        static QString select_file();

        /**
         * @brief Set the current file type for conditional visualization
         * @param type File type classification (MIT/OHNE/UNKNOWN)
         */
        void set_file_type(FileType type);

        /**
         * @brief Get the current file type
         * @return Current file type classification
         */
        [[nodiscard]] FileType get_file_type() const { return currentFileType; }

        /**
         * @brief Plot raw waveforms overlay and average pulse shapes
         *
         * Creates two side-by-side charts:
         * - Left: All waveforms overlaid with gate lines
         * - Right: Average neutron vs gamma pulse shapes (log scale)
         *   For background files, shows average waveform only.
         *
         * @param waveforms Processed waveform data
         * @param short_gate_ns Short gate duration for gate line display
         * @param total_gate_ns Total gate duration for gate line display
         */
        void plot_raw_waveforms(const std::vector<WaveformData>& waveforms, double short_gate_ns = 25.0, double total_gate_ns = 150.0);

        /**
         * @brief Plot average pulse shapes for neutron/gamma comparison
         *
         * Separates waveforms by PSD value and plots average shapes:
         * - Blue: Low PSD (gamma events)
         * - Red: High PSD (neutron events)
         *
         * For background files, plots single average waveform.
         *
         * @param waveforms Processed waveform data
         * @param psd_params Corresponding PSD parameters
         */
        void plot_average_pulse_shapes(const std::vector<WaveformData>& waveforms, const std::vector<PSDParameters>& psd_params) const;

        /**
         * @brief Plot 2D PSD scatter plot
         *
         * Creates a scatter plot of PSD vs Qtot with color-coded points.
         * Includes median PSD discrimination line.
         * Only displayed for MIT (measurement) files.
         *
         * @param psd_params PSD parameters for all events
         */
        void plot_psd_scatter(const std::vector<PSDParameters>& psd_params);

        /**
         * @brief Plot interactive PSD histogram with parameter sliders
         *
         * Creates a histogram with sliders for:
         * - Short gate duration (ns)
         * - Total gate duration (ns)
         * - Energy threshold
         * - Number of bins
         *
         * Updates in real-time as sliders are adjusted.
         * Only active for MIT (measurement) files.
         *
         * @param psd_params PSD parameters for all events
         * @param energy_threshold Initial energy threshold value
         * @param fom_result Initial FOM calculation result
         * @param short_gate_ns Initial short gate duration
         * @param total_gate_ns Initial total gate duration
         */
        void plot_psd_histogram(const std::vector<PSDParameters>& psd_params,
                                double energy_threshold = 10.0,
                                const FOMResult& fom_result = FOMResult(),
                                double short_gate_ns = 25.0,
                                double total_gate_ns = 150.0);

        /**
         * @brief Plot statistical distribution histograms
         *
         * Creates four-panel display:
         * - Amplitude distribution
         * - Rise time distribution
         * - Qtot vs Qshort scatter (MIT only)
         * - PSD distribution (MIT only)
         *
         * For background files, shows placeholder text in PSD panels.
         *
         * @param waveforms Processed waveform data
         * @param psd_params Corresponding PSD parameters
         */
        void plot_statistics(const std::vector<WaveformData>& waveforms, const std::vector<PSDParameters>& psd_params) const;

        /**
         * @brief Display analysis summary report
         *
         * Shows formatted statistics report in text widget.
         *
         * @param waveforms Processed waveform data
         * @param psd_params Corresponding PSD parameters
         */
        void show_summary(const std::vector<WaveformData>& waveforms, const std::vector<PSDParameters>& psd_params) const;

    private slots:
        /**
         * @brief Handle short gate slider change
         * @param value New slider value (nanoseconds)
         */
        void on_short_gate_changed(int value);

        /**
         * @brief Handle total gate slider change
         * @param value New slider value (nanoseconds)
         */
        void on_total_gate_changed(int value);

        /**
         * @brief Handle energy threshold slider change
         * @param value New slider value
         */
        void on_energy_threshold_changed(int value);

        /**
         * @brief Handle histogram bins slider change
         * @param value New number of bins
         */
        void on_bins_changed(int value);

        /**
         * @brief Recalculate and redraw PSD histogram
         */
        void update_psd_histogram();

        /**
         * @brief Save current plot to file
         */
        void save_plot();

        /**
         * @brief Export analysis data to CSV
         */
        void export_data();

    private:
        QTabWidget* tabWidget;

        QChartView* rawWaveformChartView;
        QChartView* avgWaveformChartView;

        QChartView* psdScatterChartView;

        QChartView* psdHistogramChartView;
        QSlider* shortGateSlider;
        QSlider* totalGateSlider;
        QSlider* energyThresholdSlider;
        QSlider* binsSlider;
        QLabel* shortGateLabel;
        QLabel* totalGateLabel;
        QLabel* energyLabel;
        QLabel* binsLabel;
        QLabel* fomLabel;

        QChartView* amplitudeHistView;
        QChartView* riseTimeHistView;
        QChartView* qtotQshortView;
        QChartView* psdDistView;

        QTextEdit* summaryText;

        std::vector<WaveformData> currentWaveforms;
        std::vector<PSDParameters> currentPsdParams;
        double currentShortGate = 25.0;
        double currentTotalGate = 150.0;
        double currentEnergyThreshold = 10.0;
        int currentBins = 100;
        FOMResult currentFomResult;
        FileType currentFileType = FileType::MIT;

        /**
         * @brief Create time axis in nanoseconds
         * @param n_samples Number of samples
         * @return Vector of time values (ns)
         */
        static std::vector<double> create_time_axis(size_t n_samples);

        /**
         * @brief Map value to blue-red color gradient
         * @param value Value to map
         * @param min_val Minimum of range
         * @param max_val Maximum of range
         * @return QColor for the value
         */
        static QColor get_color_from_value(double value, double min_val, double max_val);

        /**
         * @brief Create histogram chart from data vector
         * @param data Input values
         * @param bins Number of bins
         * @param title Chart title
         * @param xlabel X-axis label
         * @param ylabel Y-axis label
         * @param color Bar color
         * @return Created QChart
         */
        static QChart* create_value_histogram(const std::vector<double>& data,
                                              int bins,
                                              const QString& title,
                                              const QString& xlabel,
                                              const QString& ylabel,
                                              const QColor& color);

        /**
         * @brief Add vertical gate indicator lines to chart
         * @param chart Target chart
         * @param short_gate_ns Short gate position (ns)
         * @param total_gate_ns Total gate position (ns)
         * @param max_y Maximum Y value for lines
         */
        static void add_gate_lines(QChart* chart, double short_gate_ns, double total_gate_ns, double max_y);
    };
}