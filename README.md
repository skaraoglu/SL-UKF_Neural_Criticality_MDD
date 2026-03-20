<div align="center">

# Neural Criticality in Major Depressive Disorder

### Real-Time fMRI Neurofeedback & Stuart-Landau Whole-Brain Dynamics

[![R](https://img.shields.io/badge/R-≥4.2-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-≥3.9-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-Academic_Use-lightgrey)]()
[![Status](https://img.shields.io/badge/Status-Pilot_Complete-brightgreen)]()
[![Atlas](https://img.shields.io/badge/Parcellation-216_ROI-blue)]()
[![Subjects](https://img.shields.io/badge/N-19_paired-orange)]()

---

*Can we measure how far the depressed brain sits from a critical phase transition —  
and shift it back with neurofeedback?*

</div>

---

## Overview

This repository contains the analysis pipeline for a **double-blind, sham-controlled rtfMRI neurofeedback study** in unmedicated Major Depressive Disorder. The central scientific question is whether MDD resting-state brain dynamics occupy a **subcritical regime** — and whether neurofeedback training can measurably perturb that regime.

We fit a **Stuart-Landau oscillator** (the normal form of a supercritical Hopf bifurcation) to each brain region's BOLD time series via an **Unscented Kalman Filter**, estimating a per-region *bifurcation parameter* $a$ that quantifies each region's distance from the critical boundary between noise-driven and self-sustaining oscillatory dynamics.

<div align="center">

```
    a < 0              a = 0              a > 0
  ───────────────────────┼───────────────────────
  Subcritical        Critical         Supercritical
  (noise-driven)    (phase transition)  (limit cycle)
    ◂── MDD ──▸     ▲
                    │
               therapeutic
                  target
```

</div>

---

## Study Design

<table>
<tr>
<td width="50%">

**Participants**
- Unmedicated MDD (DSM-IV-TR)
- 23 enrolled → 19 paired for analysis
- 2 excluded (excessive head motion)
- 2 excluded (single-session only)

**Neurofeedback Protocol**
- Active: left amygdala upregulation
- Sham: left intraparietal sulcus (control)
- Double-blind, randomized assignment

</td>
<td width="50%">

**Acquisition**
- Siemens 3T scanner
- TR = 2.0 s, 260 volumes per session
- Resting-state fMRI pre- and post-NF
- AFNI preprocessing (motion, nuisance, bandpass 0.01–0.10 Hz)

**Parcellation**
- Primary: Schaefer-200 + Melbourne-16 subcortical (216 ROIs)
- Validation: Harvard-Oxford 110-ROI (sphere-based, independent)

</td>
</tr>
</table>

---

## Pipeline Architecture

```
┌──────────────────────────────────────────────────────────────────────────┐
│                         parcellate_219roi_v3.ipynb                       │
│  AFNI .BRIK/.HEAD  ──▸  Atlas construction  ──▸  ROI time-series CSVs    │
└───────────────────────────────────┬──────────────────────────────────────┘
                                    │
                                    ▼
┌──────────────────────────────────────────────────────────────────────────┐
│                          mdd_analysis_v3.ipynb                           │
│                                                                          │
│  ┌─────────────┐    ┌──────────────────┐    ┌────────────────────────┐   │
│  │  Stage 1    │    │  Stage 2a        │    │  Stage 2c              │   │
│  │  SL-UKF     │──▸ │  K identifiab.   │    │  mOU Lasso-MVAR        │   │
│  │  a, ω per   │    │  → PLV fallback  │    │  effective connectivity│   │
│  │  ROI        │    └──────────────────┘    └────────────────────────┘   │
│  └──────┬──────┘              │                         │                │
│         │                     ▼                         ▼                │
│         │            ┌────────────────┐       ┌──────────────────┐       │
│         │            │  PLV matrices  │       │  Brain graphs    │       │
│         │            │  (216 × 216)   │       │  Topology metrics│       │
│         │            └────────────────┘       └──────────────────┘       │
│         ▼                                                                │
│  ┌──────────────────────────────────────────────────────────────────┐    │
│  │  Hypothesis testing: H1–H4  ·  Sensitivity  ·  Cross-parcellation│    │
│  └──────────────────────────────────────────────────────────────────┘    │
└──────────────────────────────────────────────────────────────────────────┘
```

---

## The Model

The Stuart-Landau equation in complex form:

$$\dot{z} = (a + i\omega)\,z \ - \ |z|^2\,z \ + \ \eta(t)$$

Expanded to real coordinates for the UKF state-space:

$$\dot{x} = a\,x - \omega\,y - (x^2+y^2)\,x$$

$$\dot{y} = \omega\,x + a\,y - (x^2+y^2)\,y$$

where $z = x + iy$ is the analytic signal (BOLD + Hilbert transform), $a$ is the bifurcation parameter, $\omega$ is the natural frequency, and $\eta$ is complex Gaussian noise.

| Parameter | Meaning | Estimated by |
|-----------|---------|-------------|
| $a$ | Distance from critical point | UKF (Stage 1) |
| $\omega$ | Natural oscillation frequency | Hilbert instantaneous phase |
| $K$ | Inter-regional coupling | Tested → non-identifiable at TR=2s |

---

## Repository Structure

```
├── mdd_analysis_v3.ipynb          # Main analysis notebook (R kernel)
├── parcellate_219roi_v3.ipynb     # Parcellation pipeline (Python)
├── R/
│   ├── sl_models.R                # Stuart-Landau ODE definitions
│   ├── ukf_engine.R               # Unscented Kalman Filter core
│   ├── optim.R                    # Iterative & L-BFGS-B optimization
│   ├── preprocessing.R            # Smoothing & signal conditioning
│   └── constants.R                # UKF tuning constants
├── data/
│   ├── source/                    # Raw AFNI BRIK/HEAD + participants.tsv
│   └── parcellated/               # Output: ROI time-series CSVs
├── atlases/                       # Melbourne subcortical NIfTI
└── results/v3/                    # All outputs, plots, CSVs
    ├── plv/                       # PLV matrices & group comparisons
    ├── sl_stage1_results_*.csv    # Per-ROI bifurcation parameter estimates
    ├── sensitivity_v3.csv         # Sensitivity analysis variants
    └── summary_v3.csv             # Hypothesis verdicts
```

---

## Hypotheses

| | Hypothesis | Type | Test |
|---|-----------|------|------|
| **H1** | MDD resting-state dynamics are subcritical ($a < 0$) | Confirmatory | One-sample $t$-test (subject-level) |
| **H2** | Active NF shifts $a$ toward criticality vs sham | Primary clinical | Welch $t$-test on $\Delta a$ |
| **H2b** | Higher proportion of ROIs shift toward criticality | Secondary | Welch $t$-test on proportions |
| **H2c** | Network-specific $\Delta a$ effects | Exploratory | Per-network $t$-tests (Bonferroni) |
| **H2d** | NF changes spatial heterogeneity of $a$ | Exploratory | Welch $t$-test on $\Delta\text{var}(a)$ |
| **H3** | Frequency shift $\Delta\omega$ between groups | Exploratory | Welch $t$-test |
| **H4** | Network topology changes ($\sigma$, $C$, $Q$) | Exploratory | Wilcoxon rank-sum |

---

## Key Technical Features

**Stuart-Landau UKF**  — RK4-integrated sigma-point Kalman filter for joint state-parameter estimation from the BOLD analytic signal. Fixed - $\omega$ variant resolves the $a$ – $\omega$ identifiability trade-off.

**Dual-atlas validation** — Every subject-level result is independently replicated on a second atlas (Schaefer+Melbourne 216-ROI vs Harvard-Oxford 110-ROI) with no shared ROIs or processing code.

**Hilbert analytic signal** — FFT-based Hilbert transform with amplitude normalization provides both UKF observation channels ($x$, $y$) and instantaneous phase for PLV computation. Edge trimming (5 TRs) mitigates periodicity artifacts.

**Linearization coherence** — In the subcritical regime, the Stuart-Landau cubic term vanishes and the model reduces to multivariate Ornstein-Uhlenbeck. The mOU-MVAR connectivity pipeline is therefore not an independent model choice but the *linearization* of the SL dynamics confirmed by Stage 1.

---

## Requirements

<table>
<tr>
<td>

**R packages**
```
pracma, MASS, Matrix, dplyr,
tidyr, ggplot2, scales, glmnet,
igraph, parallel, zoo
```

</td>
<td>

**Python packages**
```
nibabel, nilearn, numpy,
pandas, scipy, tqdm
```

</td>
</tr>
</table>

**System:** R ≥ 4.2 · Python ≥ 3.9 · AFNI (for preprocessing only)

---

## Quick Start

```bash
# 1. Clone and set up
git clone ...
cd ...

# 2. Place source data
#    data/source/processed rest scans/    (rest1 BRIK/HEAD)
#    data/source/processed rest2 scans/   (rest2 BRIK/HEAD)
#    data/source/participants.tsv          (group assignments)
#    atlases/Tian_Subcortex_S1_3T_2009cAsym.nii.gz

# 3. Run parcellation (~20 min)
jupyter execute parcellate_219roi_v3.ipynb

# 4. Run full analysis (~5 hours)
jupyter execute mdd_analysis_v3.ipynb

# Results → results/v3/
```

---

## Citation

If you use this pipeline or build on this work, please cite:

---

<div align="center">

*Built with the [Stuart-Landau](https://en.wikipedia.org/wiki/Stuart%E2%80%93Landau_equation) normal form · [Unscented Kalman Filter](https://en.wikipedia.org/wiki/Kalman_filter#Unscented_Kalman_filter) · [Schaefer 2018](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal) + [Melbourne Subcortex](https://github.com/yetianmed/subcortex)*

</div>
