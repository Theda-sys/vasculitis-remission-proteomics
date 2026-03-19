# Towards a proteomic plasma biomarker panel for diagnosing vasculitis remission

Uwe Jerke^, Marieluise Kirchner^, Theda U.P. Bartolomaeus^,
Lovis Kling, Vojtech Kratky, Zdenka Hruskova, Vladimir Tesar,
Kai-Uwe Eckardt, Adrian Schreiber, Sofia Forslund†, Philipp Mertins†,
Ralph Kettritz†

^ Equal first authorship
† Equal senior / corresponding authorship

**Affiliations**
Experimental and Clinical Research Center (ECRC) and Department of Nephrology and Medical Intensive Care, Charité - Universitätsmedizin Berlin, Berlin, Germany; Max Delbruck Center (MDC) for Molecular Medicine, Berlin, Germany; Berlin Institute of Health (BIH) at Charité - Universitätsmedizin Berlin, Berlin, Germany; Charles University, Prague, Czech Republic; Centre for Ecological and Evolutionary Synthesis (CEES), Department of Biosciences, University of Oslo, Oslo, Norway.


---

## Overview

Active ANCA-associated vasculitis (AAV) requires intensive immunosuppressive therapy, but continued treatment after remission causes substantial morbidity. Routine biomarkers such as CRP and ANCA have limited sensitivity and specificity for identifying remission.

This repository contains the analysis code, frozen model objects, and validation outputs for a multi-stage machine learning pipeline that identifies and validates a **7-protein plasma biomarker panel**:

AHSG · CLEC3B · COMP · F9 · LRG1 · MCAM · MRC1

The panel was developed using **targeted PRM-MS proteomics** and validated in:

* an independent external cohort
* a pooled 60/40 re-split validation
* multiple internal validation procedures

---

## Repository structure

```
vasculitis-remission-proteomics/

├── data/
│   ├── splits/            # Frozen 80/20 and 60/40 train/test indices
│   ├── processed/         # Combined cohort proteomics data
│   └── feature_selection/ # LASSO and Boruta outputs

├── scripts/
│   ├── 01_feature_selection/   # MetaDeconfoundR, LASSO, Boruta
│   ├── 02_model_fitting/       # GLM, Ridge, Firth models
│   ├── 03_validation/          # Calibration, permutation, overfitting tests
│   ├── 04_subgroup/            # MPO / PR3 subgroup analyses
│   └── 05_clinical_utility/    # DCA, clinical marker comparisons, Cox models

├── models/                     # Frozen .rds model objects

├── validation_outputs/
│   ├── internal_8020/          # Berlin internal test set
│   ├── pooled_6040/            # Primary pooled validation
│   ├── concentration/          # Concentration-based model
│   ├── external_cohort/        # Berlin → Prague validation
│   ├── ridge_sensitivity/      # Ridge sensitivity analyses
│   ├── calibration/            # Calibration plots and metrics
│   └── overfitting/            # Permutation and bootstrap tests

└── supplementary_tables/
    └── extended_methods.pdf
```

---

## Reproducibility

All analysis steps are documented chronologically in:

```
supplementary_tables/extended_methods.pdf
```

Frozen model objects (`.rds`) and prediction outputs are provided so that all reported figures and tables can be reproduced without rerunning the full pipeline.

Random seeds used:

| Step          | Seed |
| ------------- | ---- |
| 80/20 split   | 123  |
| LASSO CV      | 42   |
| Boruta        | 1234 |
| 60/40 split   | 7    |
| AUC bootstrap | 7    |
| DCA bootstrap | 2023 |

---

## Data access

Patient-level proteomics data are available upon reasonable request under a data-sharing agreement in accordance with the ethical approval governing this study (Ethics board of Charité - Universitätsmedizin Berlin).

Contact:
**[tubartol@uio.no](mailto:tubartol@uio.no)**
**[ralph.kettritz@charite.de](mailto:ralph.kettritz@charite.de)**

