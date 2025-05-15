# 🧬 SARS-CoV-2 Divergence Analysis: Synonymous and Non-synonymous Rates

This R script analyzes longitudinal divergence data from chronically infected SARS-CoV-2 patients. It models the **synonymous mutation rate** (`divS`) while accounting for jumps, and compares this to patient-specific **non-synonymous rates** (`divN`). It also identifies **sweep-like events** based on sharp divergence increases.

---

## 🧠 Analysis Overview

| Part         | Description                                                                      |
|--------------|----------------------------------------------------------------------------------|
| 🔹 Part 1     | Estimate global synonymous rate and visualize per-patient trajectories.         |
| 🔹 Part 2     | Fit per-patient non-synonymous GLMs and compare to synonymous rate (divN/divS). |
| 🔹 Part 3     | Detect elevated divergence jumps indicative of selective sweeps.                |

## Structure

1. **Part 1: Synonymous GLM**
   - Load and process synonymous divergence data.
   - Fit a GLM with an identity link.
   - Bootstrap confidence intervals for the global synonymous rate.
   - Optional: Repeat GLM excluding late timepoints.

2. **Part 2: Non-synonymous GLM**
   - Loop over patient files.
   - Fit individual GLMs for non-synonymous divergence.
   - Compare slopes to global synonymous rate (divN/divS calculation).
   - Output a plot showing per-patient rates.

3. **Part 3: Sweep Detection**
   - Identify divergence “jumps” exceeding 2× divS.
   - Flag potential sweep-like or demographic shift events.
   - Plot per-patient divergence colored by jump status.
---

## 🔧 Required R Libraries

```r
library(dplyr)
library(tidyverse)
library(broom)
library(forcats)
library(car)
```

---

## 📂 Input Files

| File/Folder               | Description                                                 |
|--------------------------|-------------------------------------------------------------|
| `synonymous_data.csv`     | Longitudinal synonymous divergence data across patients.    |
| `non_syn_data/`           | Folder of CSVs with per-patient non-synonymous data.        |

---

## ⚙️ Key Parameters

| Parameter       | Purpose                                                               |
|----------------|-----------------------------------------------------------------------|
| `rate_thresh`   | Threshold to define divergence jumps (default: 1e-5).                |
| `max_days`      | Optional time filter for GLM fitting (e.g., `<=180` days).           |
| `jump_threshold`| Multiplier over `divS` to flag unusually fast `divN` jumps (e.g., 2).|

---

## 📈 Output Summary

- 🧮 **Global synonymous mutation rate** with confidence intervals (Wald + bootstrap).
- 🧪 **Per-patient divN/divS estimates** and significance markers.
- 🚩 **Sweep-like event detection** via elevated divergence rates.
- 📊 **Plots**:
  - Scatterplots of divergence trajectories.
  - Per-patient divN/divS plot with confidence intervals.
  - Sweep detection colored by rate elevation.

---

## 🧪 Reproducibility

- Fixed random seed via `set.seed(123)` ensures bootstrap consistency.
