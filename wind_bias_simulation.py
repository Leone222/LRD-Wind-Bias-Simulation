"""
Saveliev Correction: Kinematic Bias Analysis Tool
=================================================

Description:
    This tool analyzes the 'Saveliev Phenomenon' — the systematic overestimation 
    of Black Hole masses due to non-gravitational wind broadening.
    
    It can run in two modes:
    1. SIMULATION: Generates synthetic JADES-like LRD data to demonstrate the effect.
    2. ANALYSIS: Accepts a user-provided CSV file to calculate the correction for real data.

Author: Alexander Saveliev
Paper: "Evidence of Systematic Non-Virial Broadening Bias in Black Hole Mass Estimates at High Redshift: Revisiting the Nature of "Little Red Dots""
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.utils import resample
import sys
import os

# --- CONFIGURATION ---
# Set to None to run Simulation. Set to 'my_data.csv' to analyze real data.
DATA_FILE = None 

# If using real data, map your column names here:
COLUMN_MAP = {
    'FWHM': 'FWHM_obs',       # Column with observed FWHM (km/s)
    'Mass_Obs': 'log_M_BH_obs', # Column with Virial Mass Estimate (log Msun)
    'Luminosity': 'log_L',    # Column with Bolometric Luminosity (log erg/s)
    'M_Star': 'log_M_star'    # Column with Stellar Mass (log Msun) - needed for Excess calc
}

def generate_synthetic_data(n_samples=50):
    """Generates JADES-like synthetic data (The Simulation)."""
    print(">>> Mode: SIMULATION (Generating synthetic LRD data)...")
    np.random.seed(100)
    
    log_M_star = np.random.uniform(8.0, 10.0, n_samples)
    log_M_BH_true = log_M_star - 3.0 + np.random.normal(0, 0.3, n_samples) 
    v_virial = 1500 * np.sqrt(10**log_M_BH_true / 10**7)
    
    # The Saveliev Effect Source
    wind_strength = np.random.exponential(scale=1000, size=n_samples) 
    fwhm_obs = np.sqrt(v_virial**2 + wind_strength**2)
    
    log_L = 44 + (log_M_BH_true - 7) + np.random.normal(0, 0.4, n_samples)
    log_M_BH_obs = 6.5 + 0.5 * (log_L - 44) + 2.0 * np.log10(fwhm_obs / 1000)
    
    expected_M_BH_from_star = log_M_star - 3.0
    mass_excess = log_M_BH_obs - expected_M_BH_from_star

    df = pd.DataFrame({
        'log_M_star': log_M_star,
        'FWHM_obs': fwhm_obs,
        'Mass_Excess': mass_excess,
        'log_L': log_L,
        'log_M_BH_obs': log_M_BH_obs
    })
    return df

def load_real_data(filepath):
    """Loads and prepares real scientific data."""
    print(f">>> Mode: REAL DATA ANALYSIS (Loading {filepath})...")
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File {filepath} not found!")
    
    raw_df = pd.read_csv(filepath)
    
    # Rename columns to match internal logic
    df = pd.DataFrame()
    try:
        df['FWHM_obs'] = raw_df[COLUMN_MAP['FWHM']]
        df['log_M_BH_obs'] = raw_df[COLUMN_MAP['Mass_Obs']]
        df['log_L'] = raw_df[COLUMN_MAP['Luminosity']]
        df['log_M_star'] = raw_df[COLUMN_MAP['M_Star']]
    except KeyError as e:
        raise KeyError(f"Column not found in CSV: {e}. Check COLUMN_MAP in script.")

    # Calculate Excess for real data
    # Assuming local relation is M_BH ~ M_star - 3.0 (Kormendy & Ho approx)
    expected_M_BH = df['log_M_star'] - 3.0
    df['Mass_Excess'] = df['log_M_BH_obs'] - expected_M_BH
    
    return df

def run_saveliev_diagnostics(df):
    """
    Core Analysis Logic. 
    This is the 'Saveliev Method' applied to ANY dataframe.
    """
    
    # --- ROBUSTNESS TESTS ---
    n_boot = 5000
    boot_corrs = []
    for _ in range(n_boot):
        sample = resample(df)
        corr, _ = stats.spearmanr(sample['FWHM_obs'], sample['Mass_Excess'])
        boot_corrs.append(corr)

    boot_ci = np.percentile(boot_corrs, [2.5, 97.5])
    orig_corr, orig_p = stats.spearmanr(df['FWHM_obs'], df['Mass_Excess'])

    # Shuffle Test
    shuffle_corrs = []
    y_shuff = df['Mass_Excess'].values.copy()
    for _ in range(n_boot):
        np.random.shuffle(y_shuff)
        corr, _ = stats.spearmanr(df['FWHM_obs'], y_shuff)
        shuffle_corrs.append(corr)

    # Partial Correlation
    slope_e, intercept_e, _, _, _ = stats.linregress(df['log_L'], df['Mass_Excess'])
    resid_excess = df['Mass_Excess'] - (slope_e * df['log_L'] + intercept_e)

    slope_f, intercept_f, _, _, _ = stats.linregress(df['log_L'], df['FWHM_obs'])
    resid_fwhm = df['FWHM_obs'] - (slope_f * df['log_L'] + intercept_f)

    partial_corr, partial_p = stats.spearmanr(resid_fwhm, resid_excess)

    # PCA
    features = ['Mass_Excess', 'FWHM_obs', 'log_L', 'log_M_star']
    X_std = (df[features] - df[features].mean()) / df[features].std()
    pca = PCA(n_components=2)
    pca.fit(X_std)
    loadings = pca.components_.T * np.sqrt(pca.explained_variance_)

    # --- VISUALIZATION ---
    fig, ax = plt.subplots(1, 2, figsize=(15, 6))

    # Plot 1: Robustness
    ax[0].hist(shuffle_corrs, bins=30, alpha=0.5, color='gray', label='Null (Shuffle)')
    ax[0].hist(boot_corrs, bins=30, alpha=0.7, color='crimson', label='Bootstrap (Data)')
    ax[0].axvline(orig_corr, color='black', linestyle='--', linewidth=2, label=f'Observed r={orig_corr:.2f}')
    ax[0].set_title(f'Robustness Check\n(Bootstrap CI: [{boot_ci[0]:.2f}, {boot_ci[1]:.2f}])')
    ax[0].set_xlabel('Spearman Correlation')
    ax[0].legend()

    # Plot 2: The "Killer" Scatter
    sc = ax[1].scatter(df['FWHM_obs'], df['Mass_Excess'], c=df['log_L'], cmap='viridis', s=80, edgecolors='k')
    sns.regplot(x='FWHM_obs', y='Mass_Excess', data=df, ax=ax[1], scatter=False, color='crimson')
    ax[1].set_xlabel('FWHM (km/s)')
    ax[1].set_ylabel('Mass Excess (dex)')
    ax[1].set_title(f'The Wind Bias: FWHM predicts Anomaly\nPartial Corr (controlling for L): r={partial_corr:.2f}, p={partial_p:.2e}')
    plt.colorbar(sc, label='log Luminosity', ax=ax[1])

    plt.tight_layout()
    plt.show()

    # --- OUTPUT ---
    print("-" * 60)
    print("       SAVELIEV CORRECTION DIAGNOSTICS")
    print("-" * 60)
    print(f"Savelyev Effect Strength (Spearman r): {orig_corr:.4f} (p={orig_p:.2e})")
    print(f"Robustness (Bootstrap 95% CI): {boot_ci}")
    print(f"Independence from Luminosity (Partial Corr): {partial_corr:.4f} (p={partial_p:.2e})")
    
    slope_sav, intercept_sav, _, _, _ = stats.linregress(np.log10(df['FWHM_obs']), df['Mass_Excess'])
    print("-" * 60)
    print("PROPOSED SAVELIEV CORRECTION FORMULA FOR THIS DATASET:")
    print(f"log M_corr = log M_obs - ({slope_sav:.2f} * log10(FWHM) + {intercept_sav:.2f})")
    print("-" * 60)
    print("PCA Loadings (Confirmation of Kinematic Origin):")
    print(pd.DataFrame(loadings[:,0], index=features, columns=['PC1_Loadings']))

if __name__ == "__main__":
    # LOGIC SWITCHER
    if DATA_FILE is None:
        # Run Simulation (Default)
        df_analysis = generate_synthetic_data()
    else:
        # Run Real Analysis
        try:
            df_analysis = load_real_data(DATA_FILE)
        except Exception as e:
            print(f"Error loading data: {e}")
            sys.exit(1)
            

    run_saveliev_diagnostics(df_analysis)
