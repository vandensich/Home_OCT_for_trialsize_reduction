# Optimizing Early Ophthalmology Clinical Trials – Code Repository

This repository contains the code and data used in the publication:

**"Optimizing Early Ophthalmology Clinical Trials: Home OCT and Modeling Can Reduce Sample Size by 20–40%"**
*Jacques Hermes, Bernhard Steiert*


---

## 📌 Overview

This repository reproduces the key results of the above publication, demonstrating how model-based approaches and home OCT monitoring can reduce sample size requirements in early ophthalmology trials. The workflow includes:

1. Model calibration
2. Visual predictive checks (VPC)
3. Virtual patient simulation
4. Bootstrap-based uncertainty estimation
5. Error propagation and final analysis

---

## 📁 Repository Structure

```
.
├── Bootstrap_fits/                  # NONMEM control file for bootstrap fitting
├── Data_plot/                       # Input patient datasets and plotting script
├── Errorpropagation_analysis/      # Bootstrap results, C models, and R scripts
├── Mathematical_model_calibration/ # Initial model calibration in NONMEM
├── Mathematical_model_vpc/         # VPC generation (NONMEM + R)
├── Virtual_patient_simulation/     # Simulation scripts for trial and SoC arms
```

---

## 🩰 Requirements

- **NONMEM** (for model calibration and bootstrapping)
- **R (>= 4.0)** 
- **C Compiler** (for model compilation in error propagation step)

---

## ⚙️ How to Use

### 1. Model Calibration

Run the NONMEM model fit using:

```bash
./Mathematical_model_calibration/Run99.ctl
```

### 2. Visual Predictive Check (VPC)

After fitting, generate VPC plots:

```bash
Rscript ./Mathematical_model_vpc/VPC_plot.r
```

### 3. Virtual Patient Simulation

These scripts simulate virtual patient outcomes using the parameter values obtained from the model fitting step (as published in the paper). The simulation results provide the patient datasets that are subsequently used in the bootstrap analysis.

Simulate outcomes under both trial and standard of care conditions:

```r
Rscript ./Virtual_patient_simulation/Simulation_StandardofCare_patients.R
Rscript ./Virtual_patient_simulation/Simulation_Trialdrug_patients.R
```

### 4. Bootstrap Analysis (Manual Preparation Required)

Use the bootstrap control file:

```bash
./Bootstrap_fits/Bootstrap_fit.ctl
```

Before running, you **must manually modify**:

- The `$DATA` field to point to the relevant patient dataset.
- The `$SIML` block to reflect the correct number of virtual patients.

> ⚠️ *This step is not fully automated.* However, all required **bootstrap results** have been **pre-generated** and are available in:

```bash
./Errorpropagation_analysis/
```

These include:

- `.csv` result and covariance matrix files
- `.RData` objects for downstream analysis

### 5. Error Propagation

To compute uncertainty and generate the final results:

```r
Rscript ./Errorpropagation_analysis/Error_prop.R
```

---

## 📊 Output

- Bootstrapped results and covariance matrices:
  ```bash
  ./Errorpropagation_analysis/*.csv
  ```
- Final uncertainty plot:
  ```bash
  ./Errorpropagation_analysis/zscore_lines.pdf
  ```

---

## 🔁 Reproducibility

This repository contains all scripts, data, and pre-generated results required to reproduce the analysis in the publication. Only minor manual edits are needed for the bootstrap stage if re-running is desired.

---

## 🧾 License

This project is licensed under the **MIT License**. See the `LICENSE` file for more information.

---

## 📚 Citation

If you use this code or reproduce results from the publication, please cite:

> Jacques Hermes, Bernhard Steiert.
> *Optimizing Early Ophthalmology Clinical Trials: Home OCT and Modeling Can Reduce Sample Size by 20–40%.*> 

---

