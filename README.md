Passive Acoustic Monitoring: Call Density Estimation

Honours Research Project (R)
Author: Selaelo Kgafela
Date: 17 December 2025

⸻

Project Overview

This repository contains the full R implementation for an Honours research project on call density estimation using passive acoustic monitoring (PAM) data. The project develops, evaluates, and compares multiple statistical strategies for estimating species call density from automatically detected acoustic signals, accounting for uncertainty, imperfect detection, and stratified sampling designs.

The methodology combines probabilistic modelling, stratified sampling, Bayesian-inspired estimation, bootstrap inference, and information-theoretic optimisation, applied to real-world bioacoustic data collected at different sampling rates (32 kHz and 48 kHz).

The code is designed to be fully reproducible, parallelised, and modular, allowing extension to additional species, binning schemes, or sampling strategies.

⸻

Key Objectives
	•	Estimate species-level call density from acoustic confidence scores
	•	Correct for missing detections and imperfect classification
	•	Compare estimation accuracy across:
	•	Multiple binning schemes
	•	Different sampling rates (32 kHz vs 48 kHz)
	•	Three distinct estimation strategies
	•	Quantify uncertainty using bootstrap and empirical variance estimates
	•	Evaluate estimator performance using relative bias and ground truth comparisons

⸻

Data Description

The analysis uses two primary categories of data:
	1.	Unvalidated (Natural State) Data
	•	Combined detection outputs from PAM systems
	•	Confidence scores transformed to logit scale
	•	Data stratified by sampling rate (32 kHz and 48 kHz)
	2.	Validated Data
	•	Subsets of detections with known outcomes
	•	Used to estimate per-bin detection probabilities
	•	Provides ground truth for bias evaluation

Missing audio clips are explicitly modelled and reintegrated to avoid bias from incomplete temporal coverage.

⸻

Methodological Pipeline

1. Data Preparation
	•	Load and clean raw PAM outputs
	•	Harmonise formats across sampling rates
	•	Transform confidence scores to logit scale
	•	Reconstruct missing temporal segments probabilistically

2. Binning and Annotation
	•	Partition logit scores into quantile-based bins
	•	Estimate bin-level positive and negative probabilities using validated data
	•	Compute bin weights for stratified aggregation

3. Density Estimation Strategies

Study-Level Estimation
	•	Baseline density estimator combining both strata
	•	Includes bootstrap-based uncertainty quantification

Strategy 1: Weighted Beta-Binomial Estimation
	•	Independent density estimates per sampling rate
	•	Uses bin-level beta smoothing and stratified bootstrapping

Strategy 2: Information-Theoretic (KL-Minimisation) Estimation
	•	Optimises site-level density via KL divergence
	•	Fully vectorised and parallelised for computational efficiency

Strategy 3: Hybrid Estimator
	•	Geometric combination of Strategy 1 and Strategy 2
	•	Designed to stabilise bias–variance trade-offs

⸻

Statistical Techniques Used
	•	Logistic (logit) transformations
	•	Quantile-based stratification
	•	Beta-binomial smoothing
	•	Stratified bootstrap resampling
	•	Monte Carlo simulation
	•	KL divergence minimisation
	•	Bias and variance decomposition
	•	Parallel computation using future and future.apply

⸻

Code Structure
	•	Data loading & preprocessing
	•	Species-specific data extraction (species() function)
	•	Quantiling and annotation (Quantile_function(), annotations())
	•	Validation subsampling (validation.data())
	•	Study-level density estimation
	•	Strategy 1, 2, and 3 implementations
	•	Bootstrap inference and uncertainty estimation
	•	Bias analysis and visualisation
	•	Automated result tables and PDF figure outputs

All major functions are modular and reusable.

⸻

Running the Analysis
	1.	Ensure all required libraries are installed
	2.	Place cleaned CSV datasets in the working directory
	3.	Adjust parallel worker settings if needed
	4.	Select:
  •	Number of bins (quantiles)
	•	Samples per bin
	•	Number of bootstrap repetitions
	5.	Run the Estimates() function
	6.	Generate plots and tables using results_output() and density.estimates()

Results are automatically saved to:
	•	Results/ (RDS outputs)
	•	figure_results/ (PDF visualisations)

⸻

Outputs
	•	Call density estimates per species
	•	Strategy-wise comparisons (32 kHz vs 48 kHz)
	•	Bootstrap and empirical standard deviations
	•	Relative bias plots
	•	Publication-ready tables and figures

⸻

Intended Audience

This codebase is intended for:
	•	Statistical ecologists
	•	Bioacoustics researchers
	•	Applied statisticians
	•	Data scientists working with imperfect detection data
	•	Researchers interested in stratified inference and uncertainty quantification

⸻

Notes
	•	The code is computationally intensive; parallel execution is strongly recommended.
	•	Random seeds are handled carefully to ensure reproducibility under parallelism.
	•	All modelling assumptions are explicitly encoded in the functions.
  •	For project data, contact me on the attached email below.

⸻

Contact

Selaelo Kgafela
arnold.l.s@outlook.com
Honours in Statistics & Data Science
University of Cape Town
