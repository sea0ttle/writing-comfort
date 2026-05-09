# What Makes a Pen Comfortable?
### A data-backed look at grip, balance, and the feel of writing

This project combines original survey data with physical pen measurements to model what drives writing comfort — and to score a catalog of fountain pens on predicted comfort fit by hand size.

The final output is an interactive Quarto/Shiny document (`writing_comfort_portfolio.qmd`) that presents the analysis and lets readers explore comfort predictions across a catalog of fountain pens.

---

## Background

Most pen comfort advice is anecdotal. This project takes a more empirical approach: recruit writers, have them rate instruments they actually own, collect physical measurements of those instruments, and model what predicts high comfort ratings.

The core finding is that **grip diameter** and **relative balance point** are the strongest physical predictors of comfort — and that the effect of balance interacts meaningfully with hand size. Fountain pens also showed a comfort advantage over other instrument types after controlling for physical dimensions.

---

## Project Structure

```
.
├── 00_run_all.R                         # Run the full pipeline end-to-end
├── 01_clean_data_and_features.R         # Data cleaning and physical feature engineering
├── 02_maxdiff_hb_and_clusters.R         # MaxDiff HB preference model and respondent clustering
├── 03_glmm_review_and_validation.R      # GLMM model fitting, comparison, and validation
├── 04_final_model_visuals.R             # Final figures and catalog prediction base
├── config.R                             # Centralized configuration (seed, Stan, GLMM settings)
├── pipeline_utils.R                     # Thin loader — sources all R/ utility modules
├── R/
│   ├── utils_io.R                       # Package management, path helpers, file I/O
│   ├── utils_data.R                     # Data cleaning, regression helpers, feature engineering
│   ├── utils_clustering.R               # PAM clustering, MaxDiff task construction, HB clustering
│   ├── utils_models.R                   # GLMM fitting, summarization, coefficient extraction
│   └── utils_catalog.R                  # Catalog prediction, binary flag helpers, cutoff selection
├── stan/
│   └── hb_maxdiff.stan                  # Hierarchical Bayes MaxDiff Stan model
├── data/
│   └── raw/                             # Raw input CSVs (not committed)
├── writing_comfort_portfolio.qmd        # Interactive Quarto/Shiny report (final output)
└── click-footnotes.html                 # HTML include for footnote behavior in the QMD
```

---

## Configuration

All tunable parameters live in `config.R` and are consumed by every numbered script:

- **`config$seed`** — global RNG seed (default: `812`)
- **`config$stan`** — Stan MCMC settings: `chains`, `iter`, `warmup`, `refresh`, `reuse_existing`
- **`config$glmm`** — GLMM optimizer: `optimizer` and `max_fun`

Change parameters here rather than in individual scripts.

---

## Pipeline Overview

The four numbered scripts run in sequence and each saves an `.rds` bundle to `analysis_pipeline_output/` for the next step to consume. Run them all at once via:

```r
source("00_run_all.R")
```

Or step through them individually in order.

### 01 — Data Cleaning and Feature Engineering

Ingests raw survey responses and a pen specification catalog. Cleans and joins them, imputes missing physical measurements hierarchically (model → product line → brand), estimates center-of-gravity and torque from proxy regression models, and clusters respondents by writing activity and instrument type using PAM clustering with Jaccard distance.

**Output:** `analysis_pipeline_output/01_clean/01_clean_data_bundle.rds`

### 02 — MaxDiff Hierarchical Bayes and Preference Clusters

Fits a Hierarchical Bayes MaxDiff model in Stan (`stan/hb_maxdiff.stan`) to estimate individual-level attribute importance shares from the best-worst scaling questions. Segments respondents into two preference clusters — *Smoothness-Focused* and *Fit-Focused* writers — based on their HB importance profiles.

**Output:** `analysis_pipeline_output/02_maxdiff/02_maxdiff_bundle.rds`

### 03 — GLMM Review and Validation

Fits and compares a family of Generalized Linear Mixed Models (binomial, respondent random intercepts) predicting high comfort ratings from physical pen features and respondent characteristics. The selected final model is:

```
high6 ~ diameter_z + hand_size * cg_rel_used_z + is_fountain + (1 | survey_id)
```

Also includes matched-sample comparisons against alternative model specifications (press force, weight, length), within/between decompositions, and a bridge analysis connecting stated MaxDiff importance to revealed behavioral evidence.

**Output:** `analysis_pipeline_output/03_review/03_glmm_review_bundle.rds`

### 04 — Final Model Visuals and Catalog Predictions

Generates the core figures (HB cluster shares, balance × hand size curves, diameter effect, fountain pen effect) and applies the final GLMM to a catalog of fountain pens to produce comfort probability scores by hand size and posting configuration.

**Output:** `analysis_pipeline_output/04_visuals/04_visual_bundle.rds` + PNG/PDF figures

---

## Data Requirements

The pipeline expects the following input files in `data/raw/` (not committed to version control):

| File | Description |
|---|---|
| `combined_data_manual_clean_5.csv` | Pen specification catalog with physical measurements |
| `jrg.csv` | Reference dataset of pens with directly measured balance points, used to train CG proxy models |
| `Pen Survey Responses - responses (6).csv` | Main survey responses including comfort ratings and respondent attributes |
| `Pen Survey Responses - comfort_items.csv` | Item-level comfort ratings linked to specific instruments |

The catalog CSV should include columns for grip diameter, body diameter, weight, cap weight, uncapped/posted lengths, material, manufacturer, model name, and product type, among others. The survey CSVs should include respondent-level covariates (hand size, press force, handwriting size, grip posture, writing session length) alongside per-instrument comfort ratings on a 1–7 scale.

`jrg.csv` is a reference dataset of pens for which center-of-gravity was physically measured in multiple configurations. It contains columns for length (not posted, closed, posted), body and grip diameter, writing weight, cap weight, and measured CG values (raw absolute and length-normalized relative) in both posted and unposted states — specifically: `LPn`, `LCp`, `LPs`, `DBa`, `D25`, `WPn`, `WCp`, `CGu`, `CGp`, `CgPn`, `CgPs`. Script 01 uses this dataset to fit four proxy regression models — linear and ridge — that predict balance point from observable dimensions. These models are then applied to the full pen catalog and to survey instruments to estimate balance where direct measurement is unavailable.

---

## Dependencies

All R package dependencies are managed via `ensure_packages()` in `R/utils_io.R`, which installs any missing packages on first run. Key packages include:

- **Data:** `tidyverse`, `readr`
- **Modeling:** `lme4`, `glmnet`, `rstan`, `performance`, `pROC`
- **Clustering:** `cluster`, `proxy`
- **Visualization:** `ggplot2`, `ggtext`, `tidytext`, `scales`
- **Report:** `quarto`, `shiny`, `gt`, `broom.mixed`, `htmltools`

For reproducible package versions, bootstrap `renv` (the `.Rprofile` activates it automatically):

```r
renv::init()   # first time only — creates renv.lock
renv::restore() # on a fresh clone
```

---

## Rendering the Report

Once the pipeline has been run and all four `.rds` bundles exist in `analysis_pipeline_output/`, render the interactive report with:

```r
quarto::quarto_render("writing_comfort_portfolio.qmd")
```

Or preview it locally:

```r
quarto::quarto_preview("writing_comfort_portfolio.qmd")
```

To deploy to shinyapps.io:

```r
library(quarto)
quarto_publish_app(server = "shinyapps.io")
```

---

## Key Findings

- **Grip diameter** was the strongest single physical predictor of comfort — wider grips were consistently associated with higher ratings across hand sizes.
- **Relative balance point** interacted with hand size: back-heavy instruments hurt comfort most for small hands, showed a flatter relationship for medium hands, and a moderate penalty for large hands.
- **Fountain pens** showed higher predicted comfort than other instrument types at equivalent diameter and balance, even after controlling for physical dimensions.
- **Stated vs. revealed preferences diverged:** balance ranked low in MaxDiff importance for both respondent clusters, but proved to be a strong revealed predictor in the behavioral model.
- A catalog scoring tool in the interactive report allows readers to filter fountain pens by hand size and posting configuration based on predicted comfort probability.
