# 🌌 Bayesian Global Sky Model (B-GSM) - Paper 1: Validation

<div align="center">
  <img src="bgsm_landscape_abstract.png" alt="B-GSM Landscape Abstract Visualization" width="800"/>
</div>

This repository contains all code and data for the first B-GSM validation paper:

**"The Bayesian Global Sky Model (B-GSM): Validation of a Data Driven Bayesian Simultaneous Component Separation and Calibration Algorithm for EoR Foreground Modelling"**

📄 **Paper**: [arXiv:2501.01417](https://arxiv.org/abs/2501.01417)  
👥 **Authors**: George Carter, Will Handley, Mark Ashdown, Nima Razavi-Ghods (University of Cambridge)

## 🔭 Overview

B-GSM is a novel Bayesian approach to modeling radio foregrounds at frequencies <400 MHz for Epoch of Reionization (EoR) observations. Unlike previous models that use Principal Component Analysis, B-GSM:

- 🎯 Uses **nested sampling** and **Bayesian model comparison** to determine optimal emission components
- ⚡ Performs **simultaneous component separation and calibration** 
- 📊 Provides **robust error quantification** through posterior distributions
- 🔧 Handles **calibration uncertainties** in diffuse sky surveys

## 📁 Repository Contents

### 🧮 Core Analysis Scripts

| Script | Purpose | Paper Section |
|--------|---------|---------------|
| `mk24_mixed_models_v4_mock_data_gen_comp_amps_final.py` | 🏗️ Generate synthetic dataset with calibration errors and noise | §6 "The Synthetic Dataset" |
| `mock_data_v4_mk24_mixed_models_nested_sampling_perfect_noise.py` | 🎲 Nested sampling analysis for Bayesian model comparison | §7.1 "Bayesian Evidence" |
| `pix_by_pix_mk24_final2.py` | 🧩 Pixel-by-pixel approximate marginalisation implementation | §4 "Approximate Marginalisation" |
| `gen_EDGES_beams.py` | 📡 Generate synthetic absolute temperature dataset (T vs LST curves) | §6.2 "Synthetic Absolute Temperature Dataset" |
| `test_approximation2.py` | ✅ Validate the approximation used in marginal likelihood | Appendix A |

### 📈 Visualization Scripts

| Script | Purpose |
|--------|---------|
| `plot_mean_comps_and_skys_v2.py` | 🗺️ Generate component and sky map comparisons |
| `plot_posterior_results_for_paper_mixed_model_with_crosshair_v2.py` | 📊 Create posterior analysis plots and corner plots |

### 💾 Data

- 📦 `mock_dataset_final.zip` - Complete synthetic dataset used for validation
- 📄 `2501.01417.pdf` - Downloaded paper PDF
- 📂 `2501.01417/` - Paper LaTeX source and figures

## 🏆 Key Results

The validation study demonstrates that B-GSM:

1. 🎯 **Correctly identifies model structure**: Bayesian evidence strongly favors 2-component model with curved power-law spectra
2. 📏 **Accurate parameter recovery**: Spectral parameters recovered within statistical uncertainty
3. ⚖️ **Successful calibration**: Significant reduction in RMS residuals between true and predicted T vs LST curves
4. 🗺️ **Robust sky predictions**: Posterior sky maps agree with true synthetic sky across full frequency range
5. 🧠 **Prior robustness**: Sky predictions remain stable across different prior assumptions

## 🧮 Mathematical Framework

B-GSM implements a joint likelihood function combining:
- 🌐 **Diffuse likelihood**: Spatially resolved but poorly calibrated survey maps
- 📶 **Absolute temperature likelihood**: Well-calibrated but 1D temperature vs LST measurements

Key innovation is the analytical marginalisation over component amplitudes with a computationally tractable pixel-by-pixel approximation.

## 🔧 Dependencies

The code requires:
- 🐍 Python scientific stack (numpy, scipy, matplotlib)
- 🌍 HEALPix for sky map handling
- 🎲 PolyChord for nested sampling
- 🔬 Standard astronomy libraries

## 🚀 Usage

1. 🏗️ **Generate synthetic data**: Run `mk24_mixed_models_v4_mock_data_gen_comp_amps_final.py`
2. 🧮 **Perform analysis**: Execute `mock_data_v4_mk24_mixed_models_nested_sampling_perfect_noise.py`
3. 📊 **Create plots**: Use the visualization scripts to reproduce paper figures

## 📝 Citation

If you use this code or methodology, please cite:

```bibtex
@article{Carter2025,
    title={The Bayesian Global Sky Model (B-GSM): Validation of a Data Driven Bayesian Simultaneous Component Separation and Calibration Algorithm for EoR Foreground Modelling},
    author={Carter, George and Handley, Will and Ashdown, Mark and Razavi-Ghods, Nima},
    journal={arXiv preprint arXiv:2501.01417},
    year={2025}
}
```

## 🔮 Future Work

This is the first of two papers. The second paper will present:
- 🌍 Application to real observational data
- 📦 Public release of the low-frequency sky model
- 🚀 Deployment-ready B-GSM implementation

## 📧 Contact

For questions about the methodology or code, please contact:
- 👨‍🔬 George Carter: gtc30@cam.ac.uk
- 📬 Will Handley (corresponding): [contact info]
