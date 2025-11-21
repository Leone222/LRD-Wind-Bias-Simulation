# Saveliev Correction Analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17676490.svg)](https://zenodo.org/records/17676490)

Code and statistical analysis for the paper **"Evidence of Systematic Non-Virial Broadening Bias in Black Hole Mass Estimates at High Redshift"**.

This tool demonstrates the **Saveliev Correction** — a method to disentangle wind-driven spectral line broadening from gravitational kinematics in "Little Red Dots" (LRDs).

## 🚀 Quick Start (Simulation)

To reproduce the results and figures from the paper using synthetic data:

1. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

2. **Run the simulation:**
   ```bash
   python wind_bias_simulation.py
   ```

## 📂 How to Analyze Your Own Data (CSV)

This script is designed to be applied to real observational datasets. To calculate the Savelyev Correction for your specific sample:

1. **Prepare your data:** You need a `.csv` file with columns for FWHM, Virial Mass, Luminosity, and Stellar Mass.
2. **Configure the script:** Open `wind_bias_simulation.py` and edit the **CONFIGURATION** section at the top:

   ```python
   # 1. Set your filename (change None to your file path)
   DATA_FILE = 'my_galaxy_data.csv'

   # 2. Map your column names
   COLUMN_MAP = {
       'FWHM': 'fwhm_km_s',       # Your column name for FWHM
       'Mass_Obs': 'log_mbh',     # Your column name for Virial Mass
       'Luminosity': 'log_lbol',  # Your column name for Luminosity
       'M_Star': 'log_mstar'      # Your column name for Stellar Mass
   }
   ```
3. **Run the analysis:**
   ```bash
   python wind_bias_simulation.py
   ```
   The script will output the specific correction formula and diagnostic plots for your data.

## 📄 Citation

**Saveliev, A. (2025). Evidence of Systematic Non-Virial Broadening Bias in Black Hole Mass Estimates at High Redshift: Revisiting the Nature of "Little Red Dots"**

Full text and data available on Zenodo: [https://zenodo.org/records/17643994](https://zenodo.org/records/17643994)
