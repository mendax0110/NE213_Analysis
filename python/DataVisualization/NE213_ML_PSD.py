"""
NE213 ML-Based Pulse Shape Discrimination
==========================================
Neutron/gamma discrimination for NE213 liquid scintillator data
acquired with Teledyne LeCroy WaveAce 2024 (8-bit ADC, 2 GS/s).

Labeling strategy: Gaussian Mixture Model on CCM PSD vs Qtot at high energies
provides semi-supervised ground truth. ML models then extend discrimination
to low-energy regime where CCM alone fails.

Methods implemented:
    - XGBoost with ~35 handcrafted features (best accuracy/speed tradeoff)
    - 1D-CNN on raw waveforms (best low-energy performance)
    - Autoencoder for noise/pileup anomaly rejection

References:
    - Binda (2016), Uppsala thesis diva2:925705, Eq.8 CCM, Eq.9 FOM
    - Baselga et al. (2024), Appl.Sci. 14(13):5532
    - Griffiths et al. (2020), Mach.Learn.:Sci.Tech.

Usage:
    python NE213_ML_PSD.py --data-dir /path/to/data --model ne213_model.pkl
    python NE213_ML_PSD.py --mode predict --model ne213_model.pkl --file /path/to/waveform.txt
"""

import argparse
import glob
import os
import pickle
import sys
import warnings
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy import signal
from scipy.ndimage import uniform_filter1d
from sklearn.mixture import GaussianMixture
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.metrics import (
    classification_report, confusion_matrix, roc_auc_score, roc_curve
)

warnings.filterwarnings("ignore", category=FutureWarning)

# Teledyne LeCroy WaveAce 2024: 200 MHz BW, 2 GS/s interleaved, 8-bit ADC
SAMPLE_RATE_GS = 2.0
NS_PER_SAMPLE = 1.0 / SAMPLE_RATE_GS # 0.5 ns
SAMPLES_PER_NS = SAMPLE_RATE_GS # 2.0

# NE213 scintillation time constants
NE213_FAST_TAU_NS = 3.7
NE213_SLOW_TAU_NS = 32.3

# Waveform geometry
WAVEFORM_LEN = 350 # trimmed sample count per event
BASELINE_SAMPLES = 50 # pre-trigger region for baseline estimation
PRE_TRIGGER_NS = 10 # integration starts this many ns before peak

# Default CCM gate widths (Binda recommended)
DEFAULT_SHORT_GATE_NS = 25
DEFAULT_TOTAL_GATE_NS = 90

# Multi-gate widths for feature extraction (ns)
MULTI_GATE_DELAYS = [16, 24, 32, 48, 64, 96, 128]


def ns_to_samples(ns):
    return int(ns * SAMPLES_PER_NS)


def samples_to_ns(samples):
    return samples / SAMPLES_PER_NS


#==============================================================================#
#  DATA LOADING                                                                #
#==============================================================================#
def load_waveform_file(filepath):
    """
    Load waveforms from a text file produced by the LeCroy acquisition script.
    Format: every 3rd line starts with 'CH1:' followed by comma-separated int values.
    Handles both raw 16-bit transfer values and 8-bit converted values.
    """
    waveforms = []
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    for i in range(0, len(lines), 3):
        line = lines[i].strip()
        if not line:
            continue
        if line.startswith("CH1:"):
            line = line[4:]

        tokens = [t for t in line.split(",") if t.strip()]
        tokens = tokens[:WAVEFORM_LEN]
        if len(tokens) < 20:
            continue

        raw = np.array([float(t) for t in tokens])
        waveforms.append(raw)

    print(f"  Loaded {len(waveforms)} waveforms from {os.path.basename(filepath)}")
    return waveforms


def load_directory(data_dir, pattern):
    """Load all files matching *pattern*.txt in data_dir."""
    files = sorted(glob.glob(os.path.join(data_dir, f"*{pattern}*.txt")))
    all_waveforms = []
    for fp in files:
        all_waveforms.extend(load_waveform_file(fp))
    return all_waveforms, files


def classify_file_type(filepath):
    name = os.path.basename(filepath).lower()
    if "mit" in name:
        return "mit"
    if "ohne" in name:
        return "ohne"
    return None


#============================================================================#
#  PREPROCESSING                                                             #
#============================================================================#
class Preprocessor:
    """
    Baseline correction, polarity inversion, adaptive noise rejection,
    and optional Savitzky-Golay smoothing to mitigate 8-bit quantization.
    """

    def __init__(self, sg_window=7, sg_order=3, noise_sigma=3.0):
        self.sg_window = sg_window
        self.sg_order = sg_order
        self.noise_sigma = noise_sigma

    def process(self, raw_waveforms):
        """
        Returns (clean_waveforms, rejected_indices).
        Each clean waveform is baseline-corrected, inverted, and smoothed.
        """
        # Detect whether data are 16-bit LeCroy transfer values or 8-bit....
        sample_max = np.max([np.max(np.abs(w)) for w in raw_waveforms[:100]])
        is_16bit = sample_max > 300

        clean = []
        rejected = []
        baselines = []

        # First pass is: compute median baseline noise
        stds = []
        for wf in raw_waveforms:
            bl = wf[:BASELINE_SAMPLES]
            stds.append(np.std(bl))
        median_std = np.median(stds)
        noise_thresh = median_std * self.noise_sigma

        for idx, wf in enumerate(raw_waveforms):
            if is_16bit:
                wf = wf.astype(np.float64)
            else:
                wf = wf.astype(np.float64)
                wf = np.where(wf > 128, wf - 255, wf)

            # Invert Flip from negative to positive pulses...
            wf = -wf

            # Do Baseline correction
            bl_region = wf[:BASELINE_SAMPLES]
            baseline = np.mean(bl_region)
            bl_std = np.std(bl_region)
            wf = wf - baseline

            if bl_std > noise_thresh:
                rejected.append((idx, "baseline_noise", bl_std))
                continue

            # Reject unphysical negative spikes in pre-trigger to prevent noise in the CCM shot gate
            pre = wf[:min(100, len(wf))]
            if len(pre) > 0 and np.min(pre) < -(6 * bl_std + 100):
                rejected.append((idx, "negative_spike", np.min(pre)))
                continue

            # Savitzky-Golay smoothing to reducce 8-bit quantization noise
            if len(wf) >= self.sg_window:
                wf = signal.savgol_filter(wf, self.sg_window, self.sg_order)

            # Pad or trim to waveform length
            if len(wf) < WAVEFORM_LEN:
                wf = np.pad(wf, (0, WAVEFORM_LEN - len(wf)), mode="edge")
            else:
                wf = wf[:WAVEFORM_LEN]

            clean.append(wf)
            baselines.append(baseline)

        n_total = len(raw_waveforms)
        n_rej = len(rejected)
        print(f"  Preprocessing: {n_total - n_rej}/{n_total} passed "
              f"({n_rej} rejected: "
              f"{sum(1 for r in rejected if r[1] == 'baseline_noise')} noise, "
              f"{sum(1 for r in rejected if r[1] == 'negative_spike')} spikes)")
        return np.array(clean, dtype=np.float64), rejected


#============================================================================#
#  FEATURE EXTRACTION                                                        #
#============================================================================#
class FeatureExtractor:
    """
    Extracts ~35 handcrafted features per waveform covering:
      - Multi-gate charge comparison (CCM at 7 delay widths)
      - Pulse shape: rise time, fall time, FWHM, width at 10%
      - Statistical moments: kurtosis, skewness, entropy, time centroid
      - Spectral: FFT centroid, bandwidth, rolloff
      - Wavelet energy ratios (db4, 4 levels)
      - Decay autocorrelation at multiple lags
      - Rolling window means in tail
    """

    def __init__(self):
        self.feature_names = None
        self._wavelet_available = False
        try:
            import pywt
            self._wavelet_available = True
        except ImportError:
            pass

    def extract(self, waveform):
        """Return dict of features for a single preprocessed waveform."""
        f = {}
        wf = np.asarray(waveform, dtype=np.float64)
        n = len(wf)
        peak_idx = int(np.argmax(wf))
        peak_amp = wf[peak_idx]

        f["peak_amplitude"] = peak_amp
        f["peak_index"] = peak_idx
        f["total_integral"] = float(np.sum(np.clip(wf, 0, None)))

        # Rise time (10% -> 90% of peak)
        t10 = self._threshold_crossing(wf, 0.1 * peak_amp, 0, peak_idx)
        t90 = self._threshold_crossing(wf, 0.9 * peak_amp, 0, peak_idx)
        f["rise_time"] = (t90 - t10) * NS_PER_SAMPLE if (t10 is not None and t90 is not None) else 0.0

        # Fall time (90% -> 10% on decay side)
        d90 = self._threshold_crossing_decay(wf, 0.9 * peak_amp, peak_idx)
        d10 = self._threshold_crossing_decay(wf, 0.1 * peak_amp, peak_idx)
        f["fall_time"] = (d10 - d90) * NS_PER_SAMPLE if (d90 is not None and d10 is not None) else 0.0

        # Width at half maximum
        h50_l = self._threshold_crossing(wf, 0.5 * peak_amp, 0, peak_idx)
        h50_r = self._threshold_crossing_decay(wf, 0.5 * peak_amp, peak_idx)
        f["fwhm"] = (h50_r - h50_l) * NS_PER_SAMPLE if (h50_l is not None and h50_r is not None) else 0.0

        # Width at 10% of peak
        w10_l = self._threshold_crossing(wf, 0.1 * peak_amp, 0, peak_idx)
        w10_r = self._threshold_crossing_decay(wf, 0.1 * peak_amp, peak_idx)
        f["width_10pct"] = (w10_r - w10_l) * NS_PER_SAMPLE if (w10_l is not None and w10_r is not None) else 0.0

        # Multi gate CCM (Charge comparison method)....
        pre_trigger = ns_to_samples(PRE_TRIGGER_NS)
        start = max(0, peak_idx - pre_trigger)
        for delay_ns in MULTI_GATE_DELAYS:
            gate_end = min(n, start + ns_to_samples(delay_ns))
            q_short = float(np.sum(np.clip(wf[start:gate_end], 0, None)))
            q_total = f["total_integral"]
            if q_total > 0:
                f[f"ccm_{delay_ns}ns"] = (q_total - q_short) / q_total
            else:
                f[f"ccm_{delay_ns}ns"] = 0.0

        # Stats for the tail regaion...
        tail_start = min(peak_idx + ns_to_samples(30), n)
        tail = wf[tail_start:]
        if len(tail) > 5:
            f["tail_mean"] = float(np.mean(tail))
            f["tail_std"] = float(np.std(tail))
        else:
            f["tail_mean"] = 0.0
            f["tail_std"] = 0.0

        # We'll treat positive part of waveform as a probability distribution
        wf_pos = np.clip(wf, 0, None)
        total = np.sum(wf_pos)
        if total > 0:
            p = wf_pos / total
            t_axis = np.arange(n, dtype=np.float64) * NS_PER_SAMPLE

            t_mean = float(np.sum(t_axis * p))
            f["time_centroid"] = t_mean

            t_var = float(np.sum(((t_axis - t_mean) ** 2) * p))
            t_std = np.sqrt(t_var) if t_var > 0 else 1e-9
            f["pulse_variance"] = t_var
            f["pulse_skewness"] = float(np.sum(((t_axis - t_mean) / t_std) ** 3 * p))
            f["pulse_kurtosis"] = float(np.sum(((t_axis - t_mean) / t_std) ** 4 * p))

            p_nz = p[p > 0]
            f["pulse_entropy"] = float(-np.sum(p_nz * np.log(p_nz)))
        else:
            f["time_centroid"] = 0.0
            f["pulse_variance"] = 0.0
            f["pulse_skewness"] = 0.0
            f["pulse_kurtosis"] = 0.0
            f["pulse_entropy"] = 0.0

        # FFT - Fast Fourier Transform feature on tail..
        tail_fft_region = wf[peak_idx:min(peak_idx + 100, n)]
        if len(tail_fft_region) > 10:
            fft_mag = np.abs(np.fft.rfft(tail_fft_region))
            freqs = np.fft.rfftfreq(len(tail_fft_region), d=NS_PER_SAMPLE)

            total_power = np.sum(fft_mag)
            if total_power > 0:
                f["spectral_centroid"] = float(np.sum(freqs * fft_mag) / total_power)
                sc = f["spectral_centroid"]
                f["spectral_bandwidth"] = float(
                    np.sqrt(np.sum(((freqs - sc) ** 2) * fft_mag) / total_power)
                )
                # 85% spectral rolloff
                cumsum = np.cumsum(fft_mag)
                rolloff_idx = np.searchsorted(cumsum, 0.85 * total_power)
                f["spectral_rolloff"] = float(freqs[min(rolloff_idx, len(freqs) - 1)])
            else:
                f["spectral_centroid"] = 0.0
                f["spectral_bandwidth"] = 0.0
                f["spectral_rolloff"] = 0.0
        else:
            f["spectral_centroid"] = 0.0
            f["spectral_bandwidth"] = 0.0
            f["spectral_rolloff"] = 0.0

        # Wavelet features (db4, 4 levels)
        if self._wavelet_available:
            import pywt
            try:
                coeffs = pywt.wavedec(wf, "db4", level=4)
                energies = [float(np.sum(c ** 2)) for c in coeffs]
                total_e = sum(energies) + 1e-12
                for lvl, e in enumerate(energies):
                    f[f"wavelet_e{lvl}"] = e / total_e
            except Exception:
                for lvl in range(5):
                    f[f"wavelet_e{lvl}"] = 0.0
        else:
            for lvl in range(5):
                f[f"wavelet_e{lvl}"] = 0.0

        # Correlation of tail decay region...
        decay_region = wf[peak_idx:min(peak_idx + 80, n)]
        if len(decay_region) > 50:
            decay_normed = decay_region - np.mean(decay_region)
            ac_full = np.correlate(decay_normed, decay_normed, mode="full")
            ac = ac_full[len(ac_full) // 2:]
            ac = ac / (ac[0] + 1e-12)
            for lag in [5, 10, 20]:
                f[f"autocorr_lag{lag}"] = float(ac[lag]) if lag < len(ac) else 0.0
        else:
            for lag in [5, 10, 20]:
                f[f"autocorr_lag{lag}"] = 0.0

        # Max. neg. slope deay in tail
        if peak_idx + 5 < n:
            grad = np.gradient(wf[peak_idx:min(peak_idx + 60, n)])
            f["max_neg_slope"] = float(np.min(grad))
        else:
            f["max_neg_slope"] = 0.0

        return f

    def extract_batch(self, waveforms):
        """Extract features for all waveforms, return (X, feature_names)."""
        feat_dicts = [self.extract(wf) for wf in waveforms]
        self.feature_names = list(feat_dicts[0].keys())
        X = np.array([[d[k] for k in self.feature_names] for d in feat_dicts])
        X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
        return X, self.feature_names

    @staticmethod
    def _threshold_crossing(wf, threshold, start, end):
        end = min(end, len(wf))
        segment = wf[start:end]
        mask = segment >= threshold
        if not mask.any():
            return None
        return start + int(np.argmax(mask))

    @staticmethod
    def _threshold_crossing_decay(wf, threshold, peak_idx):
        segment = wf[peak_idx:]
        mask = segment <= threshold
        if not mask.any():
            return None
        return peak_idx + int(np.argmax(mask))


#============================================================================#
#  LABELING: Semi-supervised via GMM on CCM(PSD) vs Qtot                     #
#============================================================================#
def label_with_gmm(features, feature_names, energy_floor=None):
    """
    Fit a 2-component GMM in (ccm_32ns, total_integral) space on the
    high-energy subset where neutron/gamma bands are separable.
    Returns labels: 0=gamma, 1=neutron, -1=ambiguous (low confidence).

    The CCM feature at ~32ns is near-optimal for NE213 per Binda.
    Using GMM avoids hardcoded PSD thresholds and adapts to your data.
    """
    ccm_idx = feature_names.index("ccm_32ns")
    qtot_idx = feature_names.index("total_integral")

    ccm = features[:, ccm_idx]
    qtot = features[:, qtot_idx]

    # label events above this charge threshold
    if energy_floor is None:
        energy_floor = np.percentile(qtot[qtot > 0], 30)

    high_e_mask = qtot > energy_floor
    if np.sum(high_e_mask) < 50:
        print("  WARNING: too few high-energy events for GMM labeling")
        return np.full(len(features), -1, dtype=int), energy_floor

    X_gmm = np.column_stack([ccm[high_e_mask], np.log1p(qtot[high_e_mask])])
    X_gmm_scaled = StandardScaler().fit_transform(X_gmm)

    gmm = GaussianMixture(n_components=2, covariance_type="full", random_state=42)
    gmm.fit(X_gmm_scaled)
    probs = gmm.predict_proba(X_gmm_scaled)
    cluster_ids = gmm.predict(X_gmm_scaled)

    # Assign gamma=lower-ccm cluster, neutron=higher-ccm cluster
    cluster_means_ccm = [np.mean(ccm[high_e_mask][cluster_ids == c]) for c in range(2)]
    gamma_cluster = int(np.argmin(cluster_means_ccm))
    neutron_cluster = 1 - gamma_cluster

    labels = np.full(len(features), -1, dtype=int)
    j = 0
    for i in range(len(features)):
        if high_e_mask[i]:
            confidence = np.max(probs[j])
            if confidence > 0.85:
                if cluster_ids[j] == gamma_cluster:
                    labels[i] = 0
                else:
                    labels[i] = 1
            j += 1

    n_gamma = np.sum(labels == 0)
    n_neutron = np.sum(labels == 1)
    n_ambig = np.sum(labels == -1)
    print(f"  GMM labeling: {n_gamma} gamma, {n_neutron} neutron, "
          f"{n_ambig} ambiguous/low-E (energy_floor={energy_floor:.0f})")
    return labels, energy_floor


def label_background(n_samples):
    """Background (ohne) events get label=2."""
    return np.full(n_samples, 2, dtype=int)


#============================================================================#
#  PILEUP / ANOMALY REJECTION via Autoencoder                                #
#============================================================================#
def build_autoencoder(input_len):
    """
    1D convolutional autoencoder for anomaly detection.
    Trained on clean pulses; high reconstruction error flags pileup/noise.
    Returns None if tensorflow is unavailable.

    Pads input to nearest multiple of 4 so MaxPool(2)x2 and
    Conv1DTranspose(stride=2)x2 round-trip without dimension loss.
    """
    try:
        import tensorflow as tf
        from tensorflow import keras
        from tensorflow.keras import layers
    except ImportError:
        return None

    import math
    padded_len = math.ceil(input_len / 4) * 4
    pad_right = padded_len - input_len

    inp = keras.Input(shape=(input_len, 1))

    # Pad to multiple of 4
    x = layers.ZeroPadding1D(padding=(0, pad_right))(inp) if pad_right > 0 else inp

    # Encoder
    x = layers.Conv1D(32, 7, activation="relu", padding="same")(x)
    x = layers.MaxPooling1D(2)(x)
    x = layers.Conv1D(16, 5, activation="relu", padding="same")(x)
    x = layers.MaxPooling1D(2)(x)

    # Decoder
    x = layers.Conv1DTranspose(16, 5, strides=2, activation="relu", padding="same")(x)
    x = layers.Conv1DTranspose(32, 7, strides=2, activation="relu", padding="same")(x)
    x = layers.Conv1D(1, 3, activation="linear", padding="same")(x)

    # Crop back to original length
    if pad_right > 0:
        x = layers.Cropping1D(cropping=(0, pad_right))(x)

    autoencoder = keras.Model(inp, x, name="pileup_ae")
    autoencoder.compile(optimizer="adam", loss="mse")
    return autoencoder


def train_autoencoder(autoencoder, waveforms, epochs=30, batch_size=64):
    """Train autoencoder on preprocessed waveforms. Returns anomaly threshold."""
    if autoencoder is None:
        return None

    X = waveforms[:, :, np.newaxis] if waveforms.ndim == 2 else waveforms
    # Normalize per-sample
    X = X / (np.max(np.abs(X), axis=1, keepdims=True) + 1e-12)

    autoencoder.fit(X, X, epochs=epochs, batch_size=batch_size, validation_split=0.1, verbose=0)

    recon = autoencoder.predict(X, verbose=0)
    mse = np.mean((X - recon) ** 2, axis=(1, 2))
    # Threshold: 95th percentile of reconstruction error
    threshold = np.percentile(mse, 95)
    return threshold


def flag_anomalies(autoencoder, waveforms, threshold):
    """Return boolean mask: True = clean, False = anomaly (pileup/noise)."""
    if autoencoder is None or threshold is None:
        return np.ones(len(waveforms), dtype=bool)

    X = waveforms[:, :, np.newaxis] if waveforms.ndim == 2 else waveforms
    X = X / (np.max(np.abs(X), axis=1, keepdims=True) + 1e-12)
    recon = autoencoder.predict(X, verbose=0)
    mse = np.mean((X - recon) ** 2, axis=(1, 2))
    return mse <= threshold


#============================================================================#
#  ML MODELS                                                                 #
#============================================================================#
def train_xgboost(X_train, y_train, X_val, y_val):
    """
    XGBoost classifier: best accuracy/speed tradeoff for PSD.
    Falls back to sklearn GradientBoosting if xgboost is unavailable.
    """
    try:
        from xgboost import XGBClassifier
        model = XGBClassifier(
            n_estimators=500,
            max_depth=6,
            learning_rate=0.05,
            subsample=0.8,
            colsample_bytree=0.8,
            eval_metric="mlogloss",
            random_state=42,
            use_label_encoder=False,
            verbosity=0,
        )
        model.fit(
            X_train, y_train,
            eval_set=[(X_val, y_val)],
            verbose=False,
        )
    except ImportError:
        from sklearn.ensemble import GradientBoostingClassifier
        print("  xgboost not installed, using sklearn GradientBoosting")
        model = GradientBoostingClassifier(
            n_estimators=300,
            max_depth=5,
            learning_rate=0.05,
            subsample=0.8,
            random_state=42,
        )
        model.fit(X_train, y_train)
    return model


def build_cnn(input_len, n_classes=3):
    """
    1D-CNN for raw-waveform classification.
    Architecture from Griffiths et al. (2020) adapted for 350-sample NE213 pulses.
    Returns None if tensorflow is unavailable.
    """
    try:
        from tensorflow import keras
        from tensorflow.keras import layers
    except ImportError:
        return None

    model = keras.Sequential([
        layers.Input(shape=(input_len, 1)),
        layers.Conv1D(32, 7, padding="same"),
        layers.BatchNormalization(),
        layers.ReLU(),
        layers.Conv1D(64, 5, padding="same"),
        layers.BatchNormalization(),
        layers.ReLU(),
        layers.MaxPooling1D(2),
        layers.Conv1D(128, 3, padding="same"),
        layers.BatchNormalization(),
        layers.ReLU(),
        layers.MaxPooling1D(2),
        layers.Conv1D(128, 3, padding="same"),
        layers.BatchNormalization(),
        layers.ReLU(),
        layers.GlobalAveragePooling1D(),
        layers.Dense(64, activation="relu"),
        layers.Dropout(0.3),
        layers.Dense(n_classes, activation="softmax"),
    ])
    model.compile(
        optimizer=keras.optimizers.Adam(1e-3),
        loss="sparse_categorical_crossentropy",
        metrics=["accuracy"],
    )
    return model


def train_cnn(model, X_train, y_train, X_val, y_val, epochs=60, batch_size=64):
    """Train 1D-CNN with early stopping."""
    if model is None:
        return None
    from tensorflow.keras.callbacks import EarlyStopping, ReduceLROnPlateau

    X_tr = X_train[:, :, np.newaxis]
    X_v = X_val[:, :, np.newaxis]

    # Normalize per-sample to [-1, 1]
    X_tr = X_tr / (np.max(np.abs(X_tr), axis=1, keepdims=True) + 1e-12)
    X_v = X_v / (np.max(np.abs(X_v), axis=1, keepdims=True) + 1e-12)

    callbacks = [
        EarlyStopping(patience=8, restore_best_weights=True, monitor="val_loss"),
        ReduceLROnPlateau(factor=0.5, patience=4, min_lr=1e-6),
    ]

    history = model.fit(
        X_tr, y_train,
        validation_data=(X_v, y_val),
        epochs=epochs,
        batch_size=batch_size,
        callbacks=callbacks,
        verbose=0,
    )
    return history


#============================================================================#
#  EVALUATION / METRICS                                                      #
#============================================================================#
def compute_fom(psd_values, labels):
    """
    Figure of Merit: FOM = |mu_n - mu_g| / (FWHM_n + FWHM_g)
    Binda Eq. 9. Only uses labels 0 (gamma) and 1 (neutron).
    """
    gamma_psd = psd_values[labels == 0]
    neutron_psd = psd_values[labels == 1]

    if len(gamma_psd) < 10 or len(neutron_psd) < 10:
        return 0.0

    mu_g = np.mean(gamma_psd)
    mu_n = np.mean(neutron_psd)
    # FWHM = 2.355 * sigma (Gaussian approximation)
    fwhm_g = 2.355 * np.std(gamma_psd)
    fwhm_n = 2.355 * np.std(neutron_psd)

    if (fwhm_g + fwhm_n) == 0:
        return 0.0
    return abs(mu_n - mu_g) / (fwhm_g + fwhm_n)


def evaluate_model(name, y_true, y_pred, y_proba=None, class_names=None):
    """Print classification report and return metrics dict."""
    if class_names is None:
        class_names = ["gamma", "neutron", "background"]

    present = sorted(set(y_true) | set(y_pred))
    target_names = [class_names[i] for i in present if i < len(class_names)]

    print(f"\n{'='*60}")
    print(f"  {name}")
    print(f"{'='*60}")
    print(classification_report(y_true, y_pred, target_names=target_names, digits=4))
    print("Confusion matrix:")
    print(confusion_matrix(y_true, y_pred))

    metrics = {"accuracy": np.mean(y_true == y_pred)}
    if y_proba is not None and len(present) == 2:
        try:
            metrics["auc"] = roc_auc_score(y_true, y_proba[:, 1])
        except Exception:
            pass
    return metrics


#============================================================================#
#  VISUALIZATION                                                             #
#============================================================================#
def plot_psd_vs_energy(features, labels, feature_names, title="PSD vs Energy"):
    """2D scatter: CCM PSD (32ns) vs total integral, colored by label."""
    ccm_idx = feature_names.index("ccm_32ns")
    qtot_idx = feature_names.index("total_integral")

    fig, ax = plt.subplots(figsize=(10, 7))
    colors = {0: "blue", 1: "red", 2: "gray", -1: "lightgray"}
    names = {0: "gamma", 1: "neutron", 2: "background", -1: "unlabeled"}

    for lbl in sorted(set(labels)):
        mask = labels == lbl
        ax.scatter(
            features[mask, qtot_idx],
            features[mask, ccm_idx],
            c=colors.get(lbl, "black"),
            label=names.get(lbl, str(lbl)),
            alpha=0.4, s=8, edgecolors="none",
        )

    ax.set_xlabel("Qtot (total charge / energy proxy)")
    ax.set_ylabel("PSD = CCM @ 32 ns")
    ax.set_title(title)
    ax.legend(markerscale=3)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    return fig


def plot_psd_histogram(psd_values, labels, title="PSD Distribution"):
    """Histogram of PSD values split by class."""
    fig, ax = plt.subplots(figsize=(9, 5))
    for lbl, name, color in [(0, "gamma", "blue"), (1, "neutron", "red")]:
        vals = psd_values[labels == lbl]
        if len(vals) > 0:
            ax.hist(vals, bins=80, alpha=0.6, label=name, color=color, edgecolor="black", linewidth=0.3)

    fom = compute_fom(psd_values, labels)
    ax.set_xlabel("PSD (CCM @ 32ns)")
    ax.set_ylabel("Counts")
    ax.set_title(f"{title}  |  FOM = {fom:.3f}")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    return fig


def plot_average_pulses(waveforms, labels, title="Average Pulse Shapes"):
    """Average waveforms per class on log scale."""
    fig, ax = plt.subplots(figsize=(10, 5))
    t = np.arange(WAVEFORM_LEN) * NS_PER_SAMPLE

    for lbl, name, color in [(0, "gamma", "blue"), (1, "neutron", "red"), (2, "background", "gray")]:
        mask = labels == lbl
        if np.sum(mask) > 0:
            avg = np.mean(waveforms[mask], axis=0)
            ax.plot(t, avg, color=color, label=f"{name} (n={np.sum(mask)})", linewidth=1.5)

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Amplitude")
    ax.set_title(title)
    ax.set_yscale("symlog", linthresh=1)
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    return fig


def plot_feature_importance(model, feature_names, top_n=15):
    """Bar chart of top feature importances from tree model."""
    try:
        imp = model.feature_importances_
    except AttributeError:
        return None

    idx = np.argsort(imp)[::-1][:top_n]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.barh(range(len(idx)), imp[idx][::-1], color="steelblue")
    ax.set_yticks(range(len(idx)))
    ax.set_yticklabels([feature_names[i] for i in idx][::-1])
    ax.set_xlabel("Feature Importance")
    ax.set_title("Top Feature Importances (XGBoost)")
    plt.tight_layout()
    return fig


def plot_confusion_matrices(y_true, preds_dict, class_names=None):
    """Side-by-side confusion matrices for multiple models."""
    if class_names is None:
        class_names = ["gamma", "neutron", "background"]

    n = len(preds_dict)
    fig, axes = plt.subplots(1, n, figsize=(6 * n, 5))
    if n == 1:
        axes = [axes]

    for ax, (name, y_pred) in zip(axes, preds_dict.items()):
        cm = confusion_matrix(y_true, y_pred)
        im = ax.imshow(cm, interpolation="nearest", cmap="Blues")
        ax.set_title(name)
        present = sorted(set(y_true) | set(y_pred))
        tick_labels = [class_names[i] for i in present if i < len(class_names)]
        ax.set_xticks(range(len(tick_labels)))
        ax.set_xticklabels(tick_labels, rotation=45)
        ax.set_yticks(range(len(tick_labels)))
        ax.set_yticklabels(tick_labels)
        ax.set_xlabel("Predicted")
        ax.set_ylabel("True")

        for i in range(cm.shape[0]):
            for j in range(cm.shape[1]):
                ax.text(j, i, str(cm[i, j]), ha="center", va="center",
                        color="white" if cm[i, j] > cm.max() / 2 else "black")

    plt.tight_layout()
    return fig


def plot_roc_curves(y_true, probas_dict, class_names=None):
    """ROC curves for binary neutron/gamma (ignore background)."""
    if class_names is None:
        class_names = ["gamma", "neutron"]

    # Filter to only gamma/neutron
    mask = (y_true == 0) | (y_true == 1)
    if np.sum(mask) < 10:
        return None

    fig, ax = plt.subplots(figsize=(7, 6))
    for name, proba in probas_dict.items():
        if proba is None:
            continue
        p = proba[mask]
        yt = y_true[mask]
        if p.ndim == 2 and p.shape[1] >= 2:
            score = p[:, 1]
        else:
            continue
        try:
            fpr, tpr, _ = roc_curve(yt, score)
            auc_val = roc_auc_score(yt, score)
            ax.plot(fpr, tpr, label=f"{name} (AUC={auc_val:.4f})")
        except Exception:
            pass

    ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("ROC: Neutron vs Gamma")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    return fig


#============================================================================#
#  PIPELINE                                                                  #
#============================================================================#
class NE213MLPipeline:
    """
    End-to-end pipeline: load -> preprocess -> label -> train -> evaluate.
    Supports training mode (requires mit+ohne data) and prediction mode
    (loads a saved model and classifies new data).
    """

    def __init__(self, data_dir=None, short_gate_ns=DEFAULT_SHORT_GATE_NS,
                 total_gate_ns=DEFAULT_TOTAL_GATE_NS):
        self.data_dir = data_dir
        self.short_gate_ns = short_gate_ns
        self.total_gate_ns = total_gate_ns

        self.preprocessor = Preprocessor()
        self.extractor = FeatureExtractor()
        self.scaler = RobustScaler()

        self.xgb_model = None
        self.cnn_model = None
        self.autoencoder = None
        self.ae_threshold = None
        self.feature_names = None

        self.mit_waveforms = None
        self.ohne_waveforms = None

    # -- data loading --

    def load_data(self):
        """Load all mit (fusor-on) and ohne (background) files from data_dir."""
        print("Loading data...")
        self.mit_raw, mit_files = load_directory(self.data_dir, "mit")
        self.ohne_raw, ohne_files = load_directory(self.data_dir, "ohne")
        print(f"  Total: {len(self.mit_raw)} mit waveforms from {len(mit_files)} files")
        print(f"  Total: {len(self.ohne_raw)} ohne waveforms from {len(ohne_files)} files")

        if len(self.mit_raw) == 0:
            print("ERROR: no 'mit' files found. Need fusor-on data for discrimination.")
            sys.exit(1)

    def load_single_file(self, filepath):
        """Load a single file for prediction."""
        raw = load_waveform_file(filepath)
        waveforms, _ = self.preprocessor.process(raw)
        return waveforms

    # -- training pipeline --

    def train(self):
        """Full training pipeline: preprocess, label, extract, train models, evaluate."""
        self.load_data()

        # Preprocess
        print("\nPreprocessing mit waveforms...")
        self.mit_waveforms, _ = self.preprocessor.process(self.mit_raw)

        if len(self.ohne_raw) > 0:
            print("Preprocessing ohne waveforms...")
            self.ohne_waveforms, _ = self.preprocessor.process(self.ohne_raw)
        else:
            self.ohne_waveforms = np.array([]).reshape(0, WAVEFORM_LEN)

        # Feature extraction for mit
        print("\nExtracting features (mit)...")
        X_mit, self.feature_names = self.extractor.extract_batch(self.mit_waveforms)

        # GMM labeling on mit data
        print("\nLabeling mit waveforms via GMM...")
        labels_mit, energy_floor = label_with_gmm(X_mit, self.feature_names)

        # Autoencoder anomaly rejection (optional)
        ae = build_autoencoder(WAVEFORM_LEN)
        if ae is not None:
            print("\nTraining autoencoder for pileup detection...")
            clean_mask = labels_mit >= 0
            if np.sum(clean_mask) > 100:
                try:
                    self.ae_threshold = train_autoencoder(ae, self.mit_waveforms[clean_mask])
                    self.autoencoder = ae
                    ae_clean = flag_anomalies(ae, self.mit_waveforms, self.ae_threshold)
                    n_anom = np.sum(~ae_clean)
                    print(f"  Autoencoder flagged {n_anom} anomalies ({100*n_anom/len(ae_clean):.1f}%)")
                    # Mark anomalies as ambiguous
                    labels_mit[~ae_clean] = -1
                except Exception as e:
                    print(f"  Autoencoder training failed ({e}), skipping anomaly rejection")
                    self.autoencoder = None

        # Ohne features + labels
        if len(self.ohne_waveforms) > 0:
            print("\nExtracting features (ohne)...")
            X_ohne, _ = self.extractor.extract_batch(self.ohne_waveforms)
            labels_ohne = label_background(len(self.ohne_waveforms))
        else:
            X_ohne = np.array([]).reshape(0, X_mit.shape[1])
            labels_ohne = np.array([], dtype=int)

        # Combine and filter out ambiguous
        X_all = np.vstack([X_mit, X_ohne])
        labels_all = np.concatenate([labels_mit, labels_ohne])
        wf_all = np.vstack([self.mit_waveforms, self.ohne_waveforms]) if len(self.ohne_waveforms) > 0 else self.mit_waveforms

        # Keep only confidently labeled data for training
        labeled_mask = labels_all >= 0
        X_labeled = X_all[labeled_mask]
        y_labeled = labels_all[labeled_mask]
        wf_labeled = wf_all[labeled_mask]

        n_classes = len(set(y_labeled))
        print(f"\nTraining set: {len(X_labeled)} samples, {n_classes} classes")
        for c in sorted(set(y_labeled)):
            name = {0: "gamma", 1: "neutron", 2: "background"}.get(c, str(c))
            print(f"  {name}: {np.sum(y_labeled == c)}")

        # Train/test split
        X_train, X_test, y_train, y_test, wf_train, wf_test = train_test_split(
            X_labeled, y_labeled, wf_labeled,
            test_size=0.2, random_state=42, stratify=y_labeled,
        )

        # Scale features
        X_train_s = self.scaler.fit_transform(X_train)
        X_test_s = self.scaler.transform(X_test)

        print("\nTraining XGBoost...")
        self.xgb_model = train_xgboost(X_train_s, y_train, X_test_s, y_test)
        y_pred_xgb = self.xgb_model.predict(X_test_s)
        try:
            y_proba_xgb = self.xgb_model.predict_proba(X_test_s)
        except Exception:
            y_proba_xgb = None
        evaluate_model("XGBoost", y_test, y_pred_xgb, y_proba_xgb)

        self.cnn_model = build_cnn(WAVEFORM_LEN, n_classes=n_classes)
        y_proba_cnn = None
        if self.cnn_model is not None:
            print("\nTraining 1D-CNN...")
            train_cnn(self.cnn_model, wf_train, y_train, wf_test, y_test)

            wf_test_norm = wf_test[:, :, np.newaxis]
            wf_test_norm = wf_test_norm / (np.max(np.abs(wf_test_norm), axis=1, keepdims=True) + 1e-12)
            y_proba_cnn = self.cnn_model.predict(wf_test_norm, verbose=0)
            y_pred_cnn = np.argmax(y_proba_cnn, axis=1)
            evaluate_model("1D-CNN", y_test, y_pred_cnn, y_proba_cnn)
        else:
            y_pred_cnn = None

        # Do FOM compariosn...
        ccm_idx = self.feature_names.index("ccm_32ns")
        psd_test = X_test[:, ccm_idx]

        fom_ccm = compute_fom(psd_test, y_test)
        print(f"\n{'='*60}")
        print(f"  CCM Figure of Merit (traditional): {fom_ccm:.4f}")

        # ML-based effective FOM: use XGBoost predictions
        fom_xgb = compute_fom(psd_test, y_pred_xgb)
        print(f"  XGBoost effective FOM:              {fom_xgb:.4f}")

        if y_pred_cnn is not None:
            fom_cnn = compute_fom(psd_test, y_pred_cnn)
            print(f"  1D-CNN effective FOM:               {fom_cnn:.4f}")
        print(f"{'='*60}")
        print(f"  Reference: Binda (NE213 @ JET): FOM ~ 1.2-1.5")
        print(f"  Reference: Baselga (ML):        FOM ~ 1.04")

        # Do plots...
        figs = []
        figs.append(plot_psd_vs_energy(X_test, y_test, self.feature_names, title="Test Set: GMM Labels on PSD vs Energy"))
        figs.append(plot_psd_histogram(psd_test, y_test, title="Test Set: PSD Distribution (GMM labels)"))
        figs.append(plot_psd_vs_energy(X_test, y_pred_xgb, self.feature_names, title="Test Set: XGBoost Classification"))
        figs.append(plot_average_pulses(wf_test, y_test, title="Average Pulse Shapes by Class"))
        figs.append(plot_feature_importance(self.xgb_model, self.feature_names))

        preds = {"XGBoost": y_pred_xgb}
        probas = {"XGBoost": y_proba_xgb}
        if y_pred_cnn is not None:
            preds["1D-CNN"] = y_pred_cnn
            probas["1D-CNN"] = y_proba_cnn

        figs.append(plot_confusion_matrices(y_test, preds))
        roc_fig = plot_roc_curves(y_test, probas)
        if roc_fig is not None:
            figs.append(roc_fig)

        # Full dataset classification
        print("\nClassifying full mit dataset with XGBoost...")
        X_mit_s = self.scaler.transform(X_mit)
        full_pred = self.xgb_model.predict(X_mit_s)
        n_g = np.sum(full_pred == 0)
        n_n = np.sum(full_pred == 1)
        n_b = np.sum(full_pred == 2)
        total = len(full_pred)
        print(f"  gamma:      {n_g:6d}  ({100*n_g/total:.1f}%)")
        print(f"  neutron:    {n_n:6d}  ({100*n_n/total:.1f}%)")
        print(f"  background: {n_b:6d}  ({100*n_b/total:.1f}%)")

        figs.append(plot_psd_vs_energy(X_mit, full_pred, self.feature_names, title="Full MIT Dataset: XGBoost Classification"))

        plt.show()
        return figs

    def predict_file(self, filepath, model_path):
        """Classify events in a new file using a saved model."""
        self._load_model(model_path)
        waveforms = self.load_single_file(filepath)

        if len(waveforms) == 0:
            print("No valid waveforms after preprocessing.")
            return

        X, _ = self.extractor.extract_batch(waveforms)
        X_s = self.scaler.transform(X)

        y_pred = self.xgb_model.predict(X_s)
        try:
            y_proba = self.xgb_model.predict_proba(X_s)
        except Exception:
            y_proba = None

        total = len(y_pred)
        for lbl, name in [(0, "gamma"), (1, "neutron"), (2, "background")]:
            n = np.sum(y_pred == lbl)
            print(f"  {name}: {n} ({100*n/total:.1f}%)")

        plot_psd_vs_energy(X, y_pred, self.feature_names, title=f"Classification: {os.path.basename(filepath)}")
        plot_average_pulses(waveforms, y_pred, title=f"Average Pulses: {os.path.basename(filepath)}")
        plt.show()

    def save_model(self, path):
        """Save XGBoost model, scaler, and feature names."""
        data = {
            "scaler": self.scaler,
            "xgb_model": self.xgb_model,
            "feature_names": self.feature_names,
            "ae_threshold": self.ae_threshold,
            "short_gate_ns": self.short_gate_ns,
            "total_gate_ns": self.total_gate_ns,
        }
        with open(path, "wb") as f:
            pickle.dump(data, f)
        print(f"Model saved to {path}")

        if self.cnn_model is not None:
            cnn_path = path.replace(".pkl", "_cnn.keras")
            self.cnn_model.save(cnn_path)
            print(f"CNN model saved to {cnn_path}")

    def _load_model(self, path):
        """Load a saved model."""
        with open(path, "rb") as f:
            data = pickle.load(f)
        self.scaler = data["scaler"]
        self.xgb_model = data["xgb_model"]
        self.feature_names = data["feature_names"]
        self.ae_threshold = data.get("ae_threshold")
        self.extractor.feature_names = self.feature_names
        print(f"Model loaded from {path}")


#============================================================================#
#  CLI                                                                       #
#============================================================================#
def main():
    parser = argparse.ArgumentParser(
        description="NE213 ML-Based Pulse Shape Discrimination"
    )
    parser.add_argument("--data-dir", "--data_dir", dest="data_dir", type=str, default=".", help="Directory containing waveform .txt files")
    parser.add_argument("--mode", choices=["train", "predict"], default="train", help="'train' on mit+ohne data or 'predict' on a single file")
    parser.add_argument("--model", type=str, default="ne213_model.pkl", help="Path to save/load the model")
    parser.add_argument("--file", type=str, default=None, help="Single file to classify (predict mode)")
    parser.add_argument("--short-gate", "--short_gate", dest="short_gate", type=float, default=DEFAULT_SHORT_GATE_NS, help="Short gate width in ns (default: 25)")
    parser.add_argument("--total-gate", "--total_gate", dest="total_gate", type=float, default=DEFAULT_TOTAL_GATE_NS, help="Total gate width in ns (default: 90)")

    args = parser.parse_args()

    print("=" * 60)
    print("  NE213 ML Pulse Shape Discrimination")
    print("  Detector: NE213 liquid scintillator")
    print("  DAQ: LeCroy WaveAce 2024 (8-bit, 2 GS/s)")
    print(f"  Mode: {args.mode}")
    print("=" * 60)

    pipeline = NE213MLPipeline(
        data_dir=args.data_dir,
        short_gate_ns=args.short_gate,
        total_gate_ns=args.total_gate,
    )

    if args.mode == "train":
        pipeline.train()
        pipeline.save_model(args.model)
    elif args.mode == "predict":
        if args.file is None:
            print("ERROR: --file required in predict mode")
            sys.exit(1)
        pipeline.predict_file(args.file, args.model)


if __name__ == "__main__":
    main()
