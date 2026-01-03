import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from scipy import signal
from scipy.optimize import curve_fit
import tkinter as tk
from tkinter import filedialog
import os
import glob
from pathlib import Path
import pickle
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
import warnings
warnings.filterwarnings('ignore')

try:
    import tensorflow as tf
    from tensorflow import keras
    from tensorflow.keras import layers
    TENSORFLOW_AVAILABLE = True
except ImportError:
    TENSORFLOW_AVAILABLE = False
    print("TensorFlow not available. Using RandomForest only.")

# Teledyne Lecroy WaveAce 2024: 200 MHz, 2 GS/s
SAMPLE_RATE_GS = 2.0
NS_PER_SAMPLE = 1.0 / SAMPLE_RATE_GS
SAMPLES_PER_NS = 1.0 / NS_PER_SAMPLE

def ns_to_samples(ns):
    """Convert nanoseconds to sample indices"""
    return int(ns * SAMPLES_PER_NS)

def samples_to_ns(samples):
    """Convert sample indices to nanoseconds"""
    return samples / SAMPLES_PER_NS

class WaveformLoader:
    """Handles loading and classification of waveform files"""

    @staticmethod
    def load_waveform_file(filename):
        """
        Load waveforms from a text file
        :param filename: The path to the waveform text file
        :return: A list of waveforms, each as a list of floats
        """
        with open(filename, 'r') as f:
            lines = f.readlines()
            waveforms = []

            for i in range(0, len(lines), 3):
                line = lines[i].strip()
                if not line:
                    continue

                if line.startswith("CH1:"):
                    line = line[4:]

                data = [x for x in line.split(',') if x]
                data = data[:350]
                waveform = [(float(x) - 255 if float(x) > 128 else float(x)) for x in data]
                waveforms.append(waveform)

        print(f"Loaded {len(waveforms)} waveforms from file: {filename}")
        return waveforms

    @staticmethod
    def load_all_files(directory, pattern):
        """
        Load all files matching pattern in directory
        :param directory: The directory to search
        :param pattern: The pattern to match in filenames
        :return: A list of all waveforms and the list of files loaded
        """
        files = glob.glob(os.path.join(directory, f"*{pattern}*.txt"))
        all_waveforms = []
        for file in files:
            waveforms = WaveformLoader.load_waveform_file(file)
            all_waveforms.extend(waveforms)
        return all_waveforms, files

class PSDFeatureExtractor:
    """Extract multiple PSD features using different methods"""

    @staticmethod
    def charge_comparison_method(waveform, short_gate_ns=25, total_gate_ns=150, pre_trigger_ns=10):
        """
        Traditional Charge Comparison Method (CCM) - Binda Equation 8
        PSD = (Qtot - Qshort) / Qtot
        :param waveform: The input waveform as a list or numpy array
        :param short_gate_ns: The short gate duration in nanoseconds
        :param total_gate_ns: The total gate duration in nanoseconds
        :param pre_trigger_ns: THe pre-trigger duration in nanoseconds
        :return: A tuple (psd_cc, qshort, qtot)
        """
        max_idx = np.argmax(waveform)
        short_gate = ns_to_samples(short_gate_ns)
        total_gate = ns_to_samples(total_gate_ns)
        pre_trigger = ns_to_samples(pre_trigger_ns)

        start_idx = max(0, max_idx - pre_trigger)
        short_end = min(len(waveform), start_idx + short_gate)
        total_end = min(len(waveform), start_idx + total_gate)

        qshort = np.sum(waveform[start_idx:short_end])
        qtot = np.sum(waveform[start_idx:total_end])

        if qtot > 0:
            psd_cc = (qtot - qshort) / qtot
        else:
            psd_cc = 0

        return psd_cc, qshort, qtot

    @staticmethod
    def zero_crossing_method(waveform):
        """
        Zero Crossing Method (ZC)
        Counts zero crossings in the tail
        :param waveform: The input waveform as a list or numpy array
        :return: The zero crossing count normalized by tail length
        """
        max_idx = np.argmax(waveform)
        tail = waveform[max_idx:]

        # Normalize tail...
        if len(tail) > 1:
            tail_normalized = tail - np.mean(tail[-10:])
            zero_crossings = np.sum(np.diff(np.sign(tail_normalized)) != 0)
            return zero_crossings / len(tail)
        return 0

    @staticmethod
    def frequency_gradient_analysis(waveform):
        """
        Frequency Gradient Analysis (FGA)
        Analyzes frequency content of the decay
        :param waveform: The input waveform as a list or numpy array
        :return: A frequency gradient metric
        """
        max_idx = np.argmax(waveform)
        if max_idx + 50 < len(waveform):
            tail = waveform[max_idx:max_idx+50]

            # FFT analysis
            fft_vals = np.fft.fft(tail)
            fft_freq = np.fft.fftfreq(len(tail))

            # Power in high frequencies (indicator of fast decay)
            high_freq_power = np.sum(np.abs(fft_vals[np.abs(fft_freq) > 0.2]))
            total_power = np.sum(np.abs(fft_vals))

            if total_power > 0:
                return high_freq_power / total_power
        return 0

    @staticmethod
    def slope_decay_comparison(waveform):
        """
        Slope Decay Comparison Method (SDCC)
        Compares slopes at different decay regions
        :param waveform: The input waveform as a list or numpy array
        :return: A slope ratio metric
        """
        max_idx = np.argmax(waveform)
        max_val = waveform[max_idx]

        # Find 50% and 10% points
        threshold_50 = 0.5 * max_val
        threshold_10 = 0.1 * max_val

        idx_50 = next((i for i in range(max_idx, len(waveform)) if waveform[i] < threshold_50), None)
        idx_10 = next((i for i in range(max_idx, len(waveform)) if waveform[i] < threshold_10), None)

        if idx_50 and idx_10 and idx_10 > idx_50:
            # Early decay slope
            early_slope = (waveform[idx_50] - max_val) / (idx_50 - max_idx) if idx_50 > max_idx else 0
            # Late decay slope
            late_slope = (waveform[idx_10] - waveform[idx_50]) / (idx_10 - idx_50) if idx_10 > idx_50 else 0

            if early_slope != 0:
                return abs(late_slope / early_slope)
        return 0

    @staticmethod
    def tail_integral_ratio(waveform):
        """
        Tail Integral Ratio - ratio of late to early integral
        :param waveform: The input waveform as a list or numpy array
        :return: A tail integral ratio metric
        """
        max_idx = np.argmax(waveform)

        if max_idx + 100 < len(waveform):
            early_integral = np.sum(waveform[max_idx:max_idx+25])
            late_integral = np.sum(waveform[max_idx+25:max_idx+100])

            if early_integral > 0:
                return late_integral / early_integral
        return 0

    @staticmethod
    def extract_all_features(waveform):
        """
        Extract all PSD features for a single waveform
        :param waveform: The input waveform as a list or numpy array
        :return: A dictionary of extracted features
        """
        features = {}

        # Charge Comparison Method
        psd_cc, qshort, qtot = PSDFeatureExtractor.charge_comparison_method(waveform)
        features['psd_cc'] = psd_cc
        features['qshort'] = qshort
        features['qtot'] = qtot

        # Zero Crossing
        features['zc'] = PSDFeatureExtractor.zero_crossing_method(waveform)

        # Frequency Gradient Analysis
        features['fga'] = PSDFeatureExtractor.frequency_gradient_analysis(waveform)

        # Slope Decay Comparison
        features['sdcc'] = PSDFeatureExtractor.slope_decay_comparison(waveform)

        # Tail Integral Ratio
        features['tir'] = PSDFeatureExtractor.tail_integral_ratio(waveform)

        # Basic pulse characteristics
        features['max_amplitude'] = max(waveform)
        features['max_idx'] = np.argmax(waveform)

        # Rise time (10% to 90%)
        max_amp = max(waveform)
        max_idx = np.argmax(waveform)
        threshold_10 = 0.1 * max_amp
        threshold_90 = 0.9 * max_amp
        idx_10 = next((i for i, v in enumerate(waveform) if v > threshold_10), None)
        idx_90 = next((i for i, v in enumerate(waveform) if v > threshold_90), None)

        if idx_10 is not None and idx_90 is not None:
            features['rise_time'] = idx_90 - idx_10
        else:
            features['rise_time'] = 0

        # Decay characteristics
        if max_idx + 50 < len(waveform):
            decay_region = waveform[max_idx:max_idx+50]
            features['decay_mean'] = np.mean(decay_region)
            features['decay_std'] = np.std(decay_region)
        else:
            features['decay_mean'] = 0
            features['decay_std'] = 0

        return features

class MLClassifier:
    """Machine Learning classifier for neutron/gamma discrimination"""

    def __init__(self, use_deep_learning=True):
        """
        Initialize the ML classifier
        :param use_deep_learning: The flag to use neutral network in addition to Random Forest
        """
        self.use_deep_learning = use_deep_learning and TENSORFLOW_AVAILABLE
        self.scaler = StandardScaler()
        self.rf_model = None
        self.nn_model = None
        self.feature_names = None

    def build_neural_network(self, input_dim):
        """
        Build a neural network for classification
        :param input_dim: The number of input features
        :return: A compiled Keras model
        """
        model = keras.Sequential([
            layers.Input(shape=(input_dim,)),
            layers.Dense(64, activation='relu'),
            layers.Dropout(0.3),
            layers.Dense(32, activation='relu'),
            layers.Dropout(0.2),
            layers.Dense(16, activation='relu'),
            layers.Dense(3, activation='softmax')  # 3 classes: gamma, neutron, noise
        ])

        model.compile(
            optimizer='adam',
            loss='sparse_categorical_crossentropy',
            metrics=['accuracy']
        )

        return model

    def prepare_training_data(self, mit_waveforms, ohne_waveforms):
        """
        Prepare training data from mit and ohne files
        Labels: 0=gamma, 1=neutron, 2=background/noise
        :param mit_waveforms: The list of waveforms from 'mit' files
        :param ohne_waveforms: The list of waveforms from 'ohne' files
        :return: A tuple (X, y) of features and labels
        """
        print("\nExtracting features from waveforms...")

        # Extract features from all waveforms
        mit_features = []
        for wf in mit_waveforms:
            feat = PSDFeatureExtractor.extract_all_features(wf)
            mit_features.append(feat)

        ohne_features = []
        for wf in ohne_waveforms:
            feat = PSDFeatureExtractor.extract_all_features(wf)
            ohne_features.append(feat)

        # Convert to arrays
        self.feature_names = list(mit_features[0].keys())

        X_mit = np.array([[f[key] for key in self.feature_names] for f in mit_features])
        X_ohne = np.array([[f[key] for key in self.feature_names] for f in ohne_features])

        # Create labels
        # For 'mit' files: use PSD_CC to separate neutron (high PSD) from gamma (low PSD)
        psd_cc_values = X_mit[:, self.feature_names.index('psd_cc')]
        median_psd = np.median(psd_cc_values)

        y_mit = np.where(psd_cc_values > median_psd, 1, 0)  # 1=neutron, 0=gamma
        y_ohne = np.full(len(X_ohne), 2)  # 2=background/noise

        # Combine datasets
        X = np.vstack([X_mit, X_ohne])
        y = np.concatenate([y_mit, y_ohne])

        # Handle NaN and inf values
        X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

        return X, y

    def train(self, X, y):
        """
        Train both Random Forest and Neural Network
        :param X: The feature matrix
        :param y: The labels
        :return: A tuple (X_test, y_test) for evaluation
        """
        print("\nSplitting data for training...")
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42, stratify=y
        )

        # Scale features
        X_train_scaled = self.scaler.fit_transform(X_train)
        X_test_scaled = self.scaler.transform(X_test)

        print("\n" + "="*70)
        print("TRAINING RANDOM FOREST CLASSIFIER")
        print("="*70)

        # Train Random Forest
        self.rf_model = RandomForestClassifier(
            n_estimators=200,
            max_depth=10,
            min_samples_split=5,
            random_state=42,
            n_jobs=-1
        )
        self.rf_model.fit(X_train_scaled, y_train)

        # Evaluate Random Forest
        y_pred_rf = self.rf_model.predict(X_test_scaled)
        print("\nRandom Forest Performance:")
        print(classification_report(y_test, y_pred_rf,
                                    target_names=['Gamma', 'Neutron', 'Background']))

        print("\nConfusion Matrix:")
        print(confusion_matrix(y_test, y_pred_rf))

        # Feature importance
        feature_importance = self.rf_model.feature_importances_
        indices = np.argsort(feature_importance)[::-1]

        print("\nTop 5 Most Important Features:")
        for i in range(min(5, len(self.feature_names))):
            idx = indices[i]
            print(f"  {i+1}. {self.feature_names[idx]}: {feature_importance[idx]:.4f}")

        # Train Neural Network if available
        if self.use_deep_learning:
            print("\n" + "="*70)
            print("TRAINING NEURAL NETWORK CLASSIFIER")
            print("="*70)

            self.nn_model = self.build_neural_network(X_train_scaled.shape[1])

            history = self.nn_model.fit(
                X_train_scaled, y_train,
                validation_split=0.2,
                epochs=50,
                batch_size=32,
                verbose=0
            )

            # Evaluate Neural Network
            y_pred_nn = np.argmax(self.nn_model.predict(X_test_scaled, verbose=0), axis=1)
            print("\nNeural Network Performance:")
            print(classification_report(y_test, y_pred_nn,
                                        target_names=['Gamma', 'Neutron', 'Background']))

            print("\nConfusion Matrix:")
            print(confusion_matrix(y_test, y_pred_nn))

        return X_test_scaled, y_test

    def predict(self, waveform, use_nn=False):
        """
        Predict class for a single waveform
        :param waveform: THe input waveform as a list or numpy array
        :param use_nn: The flag to use neural network for prediction
        :return: A tuple (predicted_class, probabilities)
        """
        features = PSDFeatureExtractor.extract_all_features(waveform)
        X = np.array([[features[key] for key in self.feature_names]])
        X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
        X_scaled = self.scaler.transform(X)

        if use_nn and self.nn_model is not None:
            pred = np.argmax(self.nn_model.predict(X_scaled, verbose=0), axis=1)[0]
            proba = self.nn_model.predict(X_scaled, verbose=0)[0]
        else:
            pred = self.rf_model.predict(X_scaled)[0]
            proba = self.rf_model.predict_proba(X_scaled)[0]

        return pred, proba

    def save_model(self, filepath):
        """
        Save trained models
        :param filepath: The base filepath to save models
        :return: A tuple (predicted_class, probabilities)
        """
        model_data = {
            'scaler': self.scaler,
            'rf_model': self.rf_model,
            'feature_names': self.feature_names,
            'use_deep_learning': self.use_deep_learning
        }

        with open(filepath + '_rf.pkl', 'wb') as f:
            pickle.dump(model_data, f)

        if self.nn_model is not None:
            self.nn_model.save(filepath + '_nn.h5')

        print(f"\nModels saved to {filepath}_rf.pkl and {filepath}_nn.h5")

    def load_model(self, filepath):
        """
        Load trained models
        :param filepath: The base filepath to load models
        :return: Nothing
        """
        with open(filepath + '_rf.pkl', 'rb') as f:
            model_data = pickle.load(f)

        self.scaler = model_data['scaler']
        self.rf_model = model_data['rf_model']
        self.feature_names = model_data['feature_names']
        self.use_deep_learning = model_data['use_deep_learning']

        if self.use_deep_learning:
            try:
                self.nn_model = keras.models.load_model(filepath + '_nn.h5')
            except:
                print("Could not load neural network model")
                self.nn_model = None

        print(f"\nModels loaded from {filepath}")

class NE213MLAnalyzer:
    """Main analyzer with ML capabilities"""

    def __init__(self):
        """
        Initialize the NE213 ML Analyzer
        """
        self.mit_waveforms = []
        self.ohne_waveforms = []
        self.ml_classifier = None
        self.figures = []

    def select_and_load_files(self):
        """
        Select directory and load all mit and ohne files
        :return: A boolean indicating success
        """
        root = tk.Tk()
        root.withdraw()

        directory = filedialog.askdirectory(title="Select directory containing NE213 data files")

        if not directory:
            return False

        print("\n" + "="*70)
        print("LOADING WAVEFORM DATA")
        print("="*70)

        # Load 'mit' files (fusor on - expect neutrons)
        print("\nLoading 'mit' files (fusor ON)...")
        self.mit_waveforms, mit_files = WaveformLoader.load_all_files(directory, 'mit')
        print(f"Found {len(mit_files)} 'mit' files with {len(self.mit_waveforms)} total waveforms")

        # Load 'ohne' files (background)
        print("\nLoading 'ohne' files (background)...")
        self.ohne_waveforms, ohne_files = WaveformLoader.load_all_files(directory, 'ohne')
        print(f"Found {len(ohne_files)} 'ohne' files with {len(self.ohne_waveforms)} total waveforms")

        if len(self.mit_waveforms) == 0 or len(self.ohne_waveforms) == 0:
            print("\nERROR: Need both 'mit' and 'ohne' files for training!")
            return False

        return True

    def train_ml_models(self):
        """
        Train machine learning models
        :return: A tuple (X_test, y_test) for evaluation
        """
        print("\n" + "="*70)
        print("TRAINING MACHINE LEARNING MODELS")
        print("="*70)

        self.ml_classifier = MLClassifier(use_deep_learning=TENSORFLOW_AVAILABLE)

        X, y = self.ml_classifier.prepare_training_data(
            self.mit_waveforms,
            self.ohne_waveforms
        )

        print(f"\nTotal samples: {len(X)}")
        print(f"  Gamma-like: {np.sum(y == 0)}")
        print(f"  Neutron-like: {np.sum(y == 1)}")
        print(f"  Background: {np.sum(y == 2)}")

        X_test, y_test = self.ml_classifier.train(X, y)

        return X_test, y_test

    def analyze_mit_files_with_ml(self):
        """
        Analyze mit files using trained ML models
        :return: A tuple (classifications, probabilities)
        """
        print("\n" + "="*70)
        print("ANALYZING 'MIT' FILES WITH ML CLASSIFIER")
        print("="*70)

        classifications = []
        probabilities = []

        for wf in self.mit_waveforms:
            pred, proba = self.ml_classifier.predict(wf, use_nn=False)
            classifications.append(pred)
            probabilities.append(proba)

        classifications = np.array(classifications)
        probabilities = np.array(probabilities)

        print(f"\nClassification Results:")
        print(f"  Gamma: {np.sum(classifications == 0)} ({np.sum(classifications == 0)/len(classifications)*100:.1f}%)")
        print(f"  Neutron: {np.sum(classifications == 1)} ({np.sum(classifications == 1)/len(classifications)*100:.1f}%)")
        print(f"  Background/Noise: {np.sum(classifications == 2)} ({np.sum(classifications == 2)/len(classifications)*100:.1f}%)")

        return classifications, probabilities

    def plot_ml_results(self, classifications, probabilities):
        """
        Visualize ML classification results
        :param classifications: The array of classification results
        :param probabilities: The array of classification probabilities
        :return: Nothing
        """
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle('ML-Based Neutron/Gamma Discrimination Results', fontsize=16)

        # Plot 1: Classification distribution
        ax1 = axes[0, 0]
        labels = ['Gamma', 'Neutron', 'Background']
        counts = [np.sum(classifications == i) for i in range(3)]
        colors = ['blue', 'red', 'gray']
        ax1.bar(labels, counts, color=colors, alpha=0.7)
        ax1.set_ylabel('Count')
        ax1.set_title('Classification Distribution')
        ax1.grid(True, alpha=0.3)

        # Plot 2: Confidence distribution
        ax2 = axes[0, 1]
        max_proba = np.max(probabilities, axis=1)
        ax2.hist(max_proba, bins=50, color='green', alpha=0.7, edgecolor='black')
        ax2.set_xlabel('Confidence')
        ax2.set_ylabel('Count')
        ax2.set_title('Classification Confidence Distribution')
        ax2.grid(True, alpha=0.3)

        # Plot 3: Sample waveforms by class
        ax3 = axes[1, 0]
        time_axis = np.arange(350) * NS_PER_SAMPLE

        for class_idx, (label, color) in enumerate(zip(labels, colors)):
            class_waveforms = [self.mit_waveforms[i] for i in range(len(classifications))
                               if classifications[i] == class_idx]
            if len(class_waveforms) > 0:
                avg_wf = np.mean(class_waveforms[:min(100, len(class_waveforms))], axis=0)
                ax3.plot(time_axis[:len(avg_wf)], avg_wf, color=color, label=label, linewidth=2, alpha=0.8)

        ax3.set_xlabel('Time (ns)')
        ax3.set_ylabel('Amplitude')
        ax3.set_title('Average Waveforms by Classification')
        ax3.legend()
        ax3.grid(True, alpha=0.3)

        # Plot 4: Feature space visualization (PSD_CC vs Qtot)
        ax4 = axes[1, 1]
        features_list = [PSDFeatureExtractor.extract_all_features(wf) for wf in self.mit_waveforms]
        psd_cc_vals = [f['psd_cc'] for f in features_list]
        qtot_vals = [f['qtot'] for f in features_list]

        for class_idx, (label, color) in enumerate(zip(labels, colors)):
            mask = classifications == class_idx
            ax4.scatter(np.array(qtot_vals)[mask], np.array(psd_cc_vals)[mask],
                        c=color, label=label, alpha=0.5, s=20)

        ax4.set_xlabel('Qtot (Energy Proxy)')
        ax4.set_ylabel('PSD Parameter')
        ax4.set_title('Feature Space: Traditional PSD vs ML Classification')
        ax4.legend()
        ax4.grid(True, alpha=0.3)

        plt.tight_layout()
        self.figures.append(fig)

    def run_analysis(self):
        """
        Run complete ML-enhanced analysis
        :return: Nothing
        """
        if not self.select_and_load_files():
            print("No files loaded. Exiting.")
            return

        self.train_ml_models()
        classifications, probabilities = self.analyze_mit_files_with_ml()
        self.plot_ml_results(classifications, probabilities)

        # Save models...
        save_path = filedialog.asksaveasfilename(
            title="Save trained models as",
            defaultextension=".pkl",
            filetypes=[("Model files", "*.pkl")]
        )

        if save_path:
            # Remove extension as we'll add it, depending on model type and the usecase...
            save_path = save_path.replace('_rf.pkl', '').replace('_nn.h5', '')
            self.ml_classifier.save_model(save_path)

        plt.show()

def main():
    """Main function of the ML analyzuzer"""
    print("="*70)
    print("NE213 ML-ENHANCED NEUTRON DETECTION ANALYSIS")
    print("Farnsworth Fusion Reactor - Machine Learning Edition")
    print("="*70)
    print("\nThis tool uses machine learning to classify neutron/gamma signals")
    print("Multiple PSD methods: CCM, ZC, FGA, SDCC, TIR")
    print("Classifiers: Random Forest + Neural Network (if TensorFlow available)")
    print("\nTensorFlow available:", TENSORFLOW_AVAILABLE)
    print("="*70)

    analyzer = NE213MLAnalyzer()
    analyzer.run_analysis()

if __name__ == "__main__":
    main()