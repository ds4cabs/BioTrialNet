# BioTrialNet

**Clinical Trial Simulation & Reanalysis Framework in R**

A teaching and exploratory research toolkit for CABS data science interns to simulate, interrogate, and reanalyze oncology clinical trials using real-world design parameters.

---

## Table of Contents

1. [Project Overview](#project-overview)
2. [Background & Motivation](#background--motivation)
3. [Statistical Concepts Covered](#statistical-concepts-covered)
4. [Repository Structure](#repository-structure)
5. [Prerequisites](#prerequisites)
6. [Setup & Installation](#setup--installation)
7. [Running the Simulation](#running-the-simulation)
8. [Script Walkthrough](#script-walkthrough)
9. [Key Parameters](#key-parameters)
10. [Outputs](#outputs)
11. [Learning Exercises](#learning-exercises)
12. [Contributing Guidelines](#contributing-guidelines)
13. [References](#references)

---

## Project Overview

BioTrialNet is a hands-on R simulation framework designed for data science interns in the CABS (Computational and Applied Biostatistics) rotation. The project models clinical trial dynamics — enrollment, event accrual, censoring, and hazard ratios — to build intuition for how trial design choices and data-handling decisions shape observed results.

**What this is:**
- A pedagogical simulation engine grounded in published trial parameters
- A sandbox for sensitivity analysis and "what-if" questions
- A starting point for future interns to build their own reanalyses or simulation extensions

**What this is NOT:**
- A reconstruction of real patient-level data
- A regulatory-grade power calculation tool
- A definitive analysis of any specific trial

---

## Background & Motivation

The inaugural simulation (`harmoni2_trial_simulation.R`) is modeled after the **HARMONi-2** trial, a randomized Phase III study comparing **ivonescimab** (a bispecific PD-1 × VEGF antibody) versus **pembrolizumab** as first-line therapy in PD-L1-high non-small cell lung cancer (NSCLC).

HARMONi-2 generated substantial scientific discussion because:

- It reported a strong **PFS benefit** (HR ~0.51) with a relatively short follow-up (~15 months from first patient in)
- The **OS data were immature** (~39% OS maturity at the PFS data cut), leaving the survival benefit uncertain
- There was a notable **censoring imbalance**: approximately 19 patients in the ivonescimab arm discontinued treatment early vs. ~4 in the pembrolizumab arm — and how that discontinuation was handled analytically affects OS interpretation

This simulation lets you quantitatively explore those three forces — hazard ratio, follow-up maturity, and informative censoring — without claiming access to the true individual patient data.

**Why this matters for drug development:** Oncology asset teams, biostatisticians, and clinical scientists routinely need to assess whether a competitor's results are robust to alternative analytical assumptions. Simulation is a principled way to do this transparently.

---

## Statistical Concepts Covered

Working through this project will build familiarity with:

| Concept | Where It Appears |
|---|---|
| Exponential survival model | `rexp_from_median()`, `simulate_trial()` |
| Hazard ratio (HR) interpretation | `safe_cox_hr()`, all summaries |
| Cox proportional hazards regression | `coxph()` calls |
| Kaplan-Meier estimation | `survfit()`, `ggsurvplot()` |
| Administrative censoring | `make_followup_limit()` |
| Informative vs. non-informative censoring | `os_bad` vs. `os_proper` scenarios |
| Data maturity (event fraction) | `calc_maturity()` |
| Monte Carlo simulation / power estimation | `run_many_trials()` |
| Sensitivity analysis | Sections 9, 10, 11 |
| Piecewise enrollment | `simulate_enrollment()` |

---

## Repository Structure

```
BioTrialNet/
├── readme.md                          # This file
└── harmoni2_trial_simulation.R        # HARMONi-2 style simulation script
```

Future interns are encouraged to add new scripts following the naming convention:

```
<trial_name>_<analysis_type>.R
```

Examples:
```
keynote189_reanalysis.R
checkmate227_simulation.R
mariposa2_censoring_sensitivity.R
```

---

## Prerequisites

### R version

R >= 4.1.0 is recommended. Check your version:

```r
R.version.string
```

### Required packages

The script auto-installs missing packages on first run. The full dependency list is:

| Package | Purpose |
|---|---|
| `survival` | Kaplan-Meier and Cox regression (core survival analysis) |
| `survminer` | Publication-quality KM plots with risk tables |
| `dplyr` | Data manipulation and summarization |
| `ggplot2` | General-purpose plotting |
| `purrr` | Functional iteration for Monte Carlo loops |
| `tibble` | Modern data frames |
| `tidyr` | Pivoting simulation results for visualization |

Install manually if needed:

```r
install.packages(c("survival", "survminer", "dplyr", "ggplot2", "purrr", "tibble", "tidyr"))
```

---

## Setup & Installation

### 1. Clone the repository

```bash
git clone <repository-url>
cd BioTrialNet
```

### 2. Open in RStudio (recommended)

Open `harmoni2_trial_simulation.R` in RStudio. The script is self-contained — no additional project configuration is required.

### 3. Verify your R environment

```r
sessionInfo()
```

Confirm R >= 4.1.0 and that the packages listed above are available or will install cleanly.

---

## Running the Simulation

### Option A — Run the entire script

Source the full file to execute all 12 sections in sequence:

```r
source("harmoni2_trial_simulation.R")
```

This will:
1. Install any missing packages
2. Run one representative simulated trial
3. Print summaries to the console
4. Generate KM plots (3 plots)
5. Run 500 Monte Carlo trials and summarize average behavior
6. Run two sensitivity analyses (censoring balance, control OS assumption)
7. Search for the data-cut that yields ~39% OS maturity
8. Export three CSV datasets

**Runtime:** approximately 2–4 minutes on a standard laptop (dominated by the 500 + 300 Monte Carlo iterations).

### Option B — Run interactively, section by section

Each section is delimited by a clear comment banner (e.g., `# --- 3) Run one simulated trial ---`). Highlight and run sections individually in RStudio to inspect intermediate results.

### Option C — Modify parameters first

Open the script, edit the `simulate_trial()` call in **Section 3** with your own assumptions, then source the file. See [Key Parameters](#key-parameters) below.

---

## Script Walkthrough

### Section 0 — Package loading

Auto-detects missing packages and installs them. Safe to re-run.

### Section 1 — Helper functions

Five utility functions used throughout:

- **`simulate_enrollment(n, accrual_months)`** — draws uniform enrollment times over the accrual window, approximating a flat enrollment rate
- **`rexp_from_median(n, median_months)`** — converts a median survival assumption into an exponential rate and samples event times
- **`make_followup_limit(entry_time, cut_date_month)`** — computes each patient's maximum observable follow-up under an administrative data cut
- **`build_time_to_event(event_time, censor_time, admin_limit)`** — resolves the minimum of true event, informative censor, and administrative censor into an observed `(time, status)` pair
- **`safe_km_median(time, status)`** and **`safe_cox_hr(time, status, arm)`** — extract KM median and Cox HR from a dataset

### Section 2 — `simulate_trial()`

The core engine. Accepts all design parameters and returns a list with three data frames:

| Object | Description |
|---|---|
| `$pfs` | PFS dataset, data cut at `pfs_cut_month` |
| `$os_proper` | OS dataset assuming treatment discontinuation does NOT censor OS (correct practice) |
| `$os_bad` | OS dataset assuming treatment discontinuation DOES censor OS (illustrates a common analytical pitfall) |

### Section 3 — Single trial run

Runs `simulate_trial()` with HARMONi-2-inspired parameters. This is the "representative example" used for all plots and printouts below.

### Section 4 — Single-trial summaries

Prints a per-arm table of: N, events, censored, KM median, HR, 95% CI, and p-value for PFS, OS-proper, and OS-bad.

### Section 5 — Kaplan-Meier plots

Three KM curves with risk tables using `survminer::ggsurvplot()`.

### Section 6 — Censoring imbalance table

Directly shows how many early-discontinuation events occurred per arm and what fraction of patients they represent. This is the most direct way to see whether censoring is balanced.

### Section 7 — `run_many_trials()`

Runs `n_sim = 500` independent trials and returns a data frame of per-simulation HR estimates, p-values, and maturity fractions. Used to estimate:
- Average HR across simulations (sanity-check against your target)
- Empirical power (fraction of simulations with p < 0.05)
- Average data maturity

### Section 8 — HR distribution visualization

Histogram of PFS and OS HRs across 500 simulations, faceted by endpoint. Shows the variance around your target HR.

### Section 9 — Censoring balance sensitivity

Compares 300 simulations with imbalanced censoring (19 vs. 4) against 300 with balanced censoring (4 vs. 4). Box plots show whether the censoring asymmetry materially shifts the HR distribution.

### Section 10 — Control OS assumption sensitivity

Varies the assumed control-arm OS median across {18, 21, 24, 27} months. Shows how an overly optimistic control assumption inflates apparent power — a common planning pitfall.

### Section 11 — Maturity calibration

Grid searches over data-cut months to find the administrative cut that yields approximately 39% OS maturity, matching the reported HARMONi-2 interim.

### Section 12 — CSV export

Writes the single representative trial's three datasets to the working directory for downstream exploration in R, Python, or Excel.

---

## Key Parameters

All parameters are documented in `simulate_trial()`. The most important ones for experimentation:

| Parameter | Default | What it controls |
|---|---|---|
| `n_armA`, `n_armB` | 198, 200 | Sample size per arm |
| `pfs_cut_month` | 15 | Months from first patient in to PFS data cut |
| `os_cut_month` | 30 | Months from first patient in to OS data cut |
| `pfs_median_B` | 5.8 | Control arm PFS median (months) |
| `pfs_hr_A_vs_B` | 0.51 | True underlying PFS hazard ratio |
| `os_median_B` | 24 | Control arm OS median (months) |
| `os_hr_A_vs_B` | 0.77 | True underlying OS hazard ratio |
| `early_disc_A_n` | 19 | Number of early discontinuations in Arm A |
| `early_disc_B_n` | 4 | Number of early discontinuations in Arm B |
| `early_disc_window` | c(1, 8) | Uniform window (months) for early discontinuation times |
| `dropout_rate` | 0.015 | Background non-informative dropout rate (exponential) |
| `accrual_months` | 10 | Duration of enrollment window |

---

## Outputs

### Console output

- Simulation assumptions (derived medians)
- Per-arm summaries for PFS, OS-proper, OS-bad
- Maturity percentages
- Censoring/early-discontinuation table
- 500-trial Monte Carlo summary
- Sensitivity analysis tables
- Top 10 data-cut calibration results

### Plots (displayed in RStudio viewer / graphics device)

| Plot | Description |
|---|---|
| KM — Simulated PFS | PFS curves with risk table and log-rank p-value |
| KM — OS (proper follow-up) | OS curves assuming proper death follow-up |
| KM — OS (if discontinuation censors OS) | OS curves under the "bad" censoring scenario |
| HR distribution histograms | Faceted by PFS / OS-proper / OS-bad across 500 sims |
| Censoring sensitivity box plots | HR spread: imbalanced vs. balanced censoring |
| Control OS bar chart | Empirical p < 0.05 rate vs. assumed control OS median |

### Exported CSV files

| File | Contents |
|---|---|
| `simulated_pfs_dataset.csv` | Per-patient PFS data with arm, entry time, true event, censor, observed time, status |
| `simulated_os_proper_dataset.csv` | Per-patient OS data (proper follow-up) |
| `simulated_os_bad_dataset.csv` | Per-patient OS data (discontinuation-censored OS) |

---

## Learning Exercises

These exercises are designed to progressively deepen your understanding. Start from wherever matches your current level.

### Beginner

1. **Reproduce the baseline results.** Source the script and confirm you see HR ~0.51 for PFS and ~0.77 for OS. Understand where those numbers come from in the code.
2. **Change the follow-up window.** Set `pfs_cut_month = 20` and re-run Section 3–5. How do the KM curves change? Why?
3. **Balance the censoring.** Set `early_disc_A_n = 4` in Section 9's `scen1`. Does the PFS HR shift? Why or why not?

### Intermediate

4. **Non-exponential survival.** The current model uses exponential (constant hazard) survival. Replace `rexp_from_median()` with a Weibull sampler (`rweibull()`) and explore how a decreasing or increasing hazard changes the KM shape and the Cox HR estimate.
5. **Add a subgroup.** Extend the data frame to include a 50/50 PD-L1 high/low split. Assume the true HR is 0.45 in PD-L1 high and 0.70 in PD-L1 low. Run forest plots of subgroup HRs across 200 simulations.
6. **Simulate an interim analysis.** Add a second data cut at `pfs_cut_month = 10` (an earlier look). Use an alpha-spending function to illustrate the multiple-testing penalty.

### Advanced

7. **Reconstruct digitized KM data.** Use the `IPDfromKM` or `reconstruct` R package to digitize a published KM curve from a trial of your choice and compare it against your simulation's output.
8. **Add a post-progression treatment effect.** Many patients cross over to novel therapy after progression. Implement a simplified rank-preserving structural failure time (RPSFT) model to adjust the OS HR for crossover.
9. **Build a Shiny app.** Wrap `simulate_trial()` and `run_many_trials()` in a Shiny interface with sliders for all key parameters so non-statisticians on the team can interactively explore scenario outcomes.
10. **Bayesian update.** Incorporate a prior on the OS HR (e.g., from the HARMONi-1 study) and derive a posterior predictive distribution for OS significance at the next planned data cut.

---

## Contributing Guidelines

Each intern cohort is encouraged to leave the repository better than they found it. Follow these conventions:

### Adding a new simulation script

1. Name your script `<trial_name>_<analysis_type>.R` (all lowercase, underscores, no spaces).
2. Include a header block at the top of the file (follow the style in `harmoni2_trial_simulation.R`):
   - Trial name and phase
   - Arms and comparator
   - Primary endpoint
   - Key published results being reproduced or interrogated
   - What the simulation is designed to teach
   - `set.seed()` for reproducibility
3. Prefer self-contained scripts that install their own packages.
4. Export at least one CSV dataset so future interns can build on your simulated data.

### Code style

- Use `<-` for assignment, not `=`
- Use `snake_case` for variable and function names
- Keep functions short and single-purpose
- Add a comment only when the _why_ is non-obvious (e.g., a statistical subtlety, a workaround)
- Do not hardcode file paths; use relative paths or `here::here()`

### Documentation

- Update this README when you add a new script: add a row to the repository structure section and a brief note in References if you are citing a published trial.
- If your script answers a specific scientific question, write one sentence at the top of the script stating what the answer is (then show the work below).

### Git workflow

```bash
git checkout -b intern/<your-initials>/<brief-description>
# make your changes
git add <files>
git commit -m "add: <trial_name> simulation for <analysis_type>"
git push origin intern/<your-initials>/<brief-description>
# open a pull request for review
```

---

## References

### Trial references

- **HARMONi-2**: Zhou C, et al. *Ivonescimab combined with chemotherapy in non-small-cell lung cancer with EGFR mutations*. NEJM 2024 (HARMONi-2 primary PFS results presented at WCLC 2024 / ESMO 2024).
- **KEYNOTE-024**: Reck M, et al. *Pembrolizumab versus chemotherapy for PD-L1-positive non-small-cell lung cancer.* NEJM 2016;375:1823–1833.

### Statistical methods

- Therneau TM, Grambsch PM. *Modeling Survival Data: Extending the Cox Model.* Springer, 2000.
- Kaplan EL, Meier P. *Nonparametric estimation from incomplete observations.* JASA 1958;53:457–481.
- Schemper M, Smith TL. *A note on quantifying follow-up in studies of failure time.* Controlled Clinical Trials 1996;17:343–346.

### R packages

- Therneau T (2024). *survival: Survival Analysis*. R package. https://CRAN.R-project.org/package=survival
- Kassambara A, Kosinski M (2021). *survminer: Drawing Survival Curves using 'ggplot2'*. https://CRAN.R-project.org/package=survminer
- Wickham H, et al. (2019). *Welcome to the tidyverse.* JOSS 4(43):1686.

### Further reading for interns

- Uno H, et al. *Moving beyond the hazard ratio in quantifying the between-group difference in survival analysis.* JCO 2014;32:2380–2385. *(Why HR alone can mislead)*
- Latouche A, et al. *A commentary on the cause-specific hazard ratio versus the subdistribution hazard ratio.* Stat Med 2013. *(Competing risks context)*
- Berry SM, et al. *Bayesian Adaptive Methods for Clinical Trials.* CRC Press, 2010. *(For the Bayesian extension exercise)*

---

*BioTrialNet is maintained by the CABS Data Science Internship Program. Questions and contributions welcome — open a PR or reach out to your rotation supervisor.*
