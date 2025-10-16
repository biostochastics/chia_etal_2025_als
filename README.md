# Confounding by geography and anticoagulant compromises proposed ALS diagnostic model and biomarkers: Re-analysis of Chia et al. (2025)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-4.4.2-blue.svg)](https://www.r-project.org/)
[![Quarto](https://img.shields.io/badge/Quarto-Report-blue.svg)](reports/confounding_investigation.qmd)

> **Reanalysis of Chia et al. (2025) Olink paper on ALS prediction, demonstrating severe confounding between plasma collection tube type, geographic origin, and diagnosis in this multi-cohort ALS biomarker study.**

---

## Read the Full Analysis

**For comprehensive findings, methodology, and discussion:**
**→ [📄 Full Published Report](https://biostochastics.quarto.pub/confounding-by-geography-and-anticoagulant-compromise-proposed-als-model/)**

*(Source: [confounding_investigation.qmd](reports/confounding_investigation.qmd))*

---

## Original Study

**Chia, R., et al. (2025).** *A plasma proteomics-based candidate biomarker panel predictive of amyotrophic lateral sclerosis.* **Nature Medicine.**

**Original Data:** https://zenodo.org/records/16899551
**Data File:** `Chia_et_al_OLINK_bridge_normalized_data_merged_deidentified_updated.txt` (312 MB)

---

## Background

This repository investigates confounding bias in the Chia et al. (2025) study. The published machine learning model uses **plasma collection tube type** (HEPARIN vs EDTA) as the 2nd most important feature for ALS diagnosis. However, tube type is perfectly confounded with geographic origin (Italy/HEPARIN vs US/EDTA) and diagnosis (both in terms of class frequencies as well as the additional analysis that was based on 100% of neurological disease controls using EDTA tubes). Our analyses demonstrate that this confounding substantially inflates model performance and undermines the validity of the proposed biomarker panel.

**Key Finding:** Model performance drops 18.2% when tested with leave-country-out cross-validation, and tube type can be predicted from protein data with 99.9% accuracy. See the full report for details.

---

## Reproducing the Analysis

### Prerequisites

- **R:** Version 4.4.2 or higher
- **Quarto:** For report generation
- **System:** 16GB RAM minimum, ~2GB disk space
- **OS:** macOS, Linux, or Windows
- **Runtime:** 2-3 hours on modern hardware

### Step 1: Clone Repository

```bash
git clone https://github.com/biostochastics/chia_etal_2025_als.git
cd chia_etal_2025_als
```

### Step 2: Restore R Environment

```r
# In R console
install.packages("renv")
renv::restore()
```

This installs all required packages with exact versions from `renv.lock`.

### Step 3: Download Original Data

Download the original dataset from Chia et al. (2025):

```bash
# Create directory
mkdir -p original

# Download from Zenodo
cd original
wget https://zenodo.org/records/16899551/files/Chia_et_al_OLINK_bridge_normalized_data_merged_deidentified_updated.txt
cd ..
```

**Note:** We do not redistribute the original data to respect data sharing agreements.

### Step 4: Run Integrated Pipeline

```r
library(targets)

# View pipeline structure (optional)
tar_visnetwork()

# Execute all analyses AND generate report
tar_make()
```

**What this does:**
- Loads and processes the OLINK data
- Runs all statistical analyses (confounding quantification, ML models, differential expression)
- Generates all figures and saves them to `outputs/figures/`
- **Automatically renders the Quarto report** when analyses complete

**Output:**
- Analysis results: `_targets/objects/` and `outputs/`
- **Final report:** `reports/confounding_investigation.html`

**Note:** The report is integrated into the pipeline via `tar_quarto()` (see `_targets.R:798-805`). The report automatically rebuilds when upstream analyses change, ensuring reproducibility.

---

## Project Organization

```
R/                    # Analysis functions (data import, ML evaluation, visualizations)
reports/              # Quarto report (comprehensive findings and methodology)
outputs/              # Generated figures and tables
original/             # Original Chia et al. data (user downloads)
_targets/             # Pipeline cache and processed objects
```


---

## Citation

### This Repository

```bibtex
@software{chia_etal_2025_als.git,
  title = {Confounding by geography and anticoagulant compromises proposed ALS diagnostic model and biomarkers: Re-analysis of Chia et al. (2025)},
  author = {Kornilov, Sergey A.},
  year = {2025},
  url = {https://github.com/biostochastics/chia_etal_2025_als}
}
```

### Original Study

```bibtex
@article{chia2025plasma,
  title = {A plasma proteomics-based candidate biomarker panel predictive of amyotrophic lateral sclerosis},
  author = {Chia, Ruth and others},
  journal = {Nature Medicine},
  year = {2025},
  url = {https://zenodo.org/records/16899551}
}
```

---

## License

- **Code:** MIT License (see [LICENSE](LICENSE))
- **Documentation:** CC-BY 4.0
- **Original Data:** Subject to Chia et al. (2025) data sharing agreements

---

## Contact

**Correspondence:** sergey.kornilov@biostochastics.com
**GitHub Issues:** For technical questions about reproduction

---

*This investigation is intended as constructive methodological insight for improving ALS biomarker research, and was submitted to Nature Medicine as a commentary on a published manuscript*
