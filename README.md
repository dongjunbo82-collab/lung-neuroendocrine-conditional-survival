# Lung neuroendocrine conditional survival analysis

R source code for stage-specific conditional survival and competing mortality comparisons between pulmonary large cell neuroendocrine carcinoma (LCNEC) and small cell lung cancer (SCLC).

This repository contains statistical analysis and plotting code only. No patient-level data, derived datasets, bootstrap files, estimates, figures, manuscripts, author contact details or local execution logs are distributed.

## Requirements

The analyses were run using R 4.6.1 with `data.table`, `survival`, `splines`, `ggplot2`, `patchwork` and `svglite`. `splines` is included with R. The scripts use base R graphics as well as ggplot2; SVG export requires svglite. Package availability is checked without installing software.

Run from the repository root. SEER inputs require independent authorized access. Set environment variables `SEER_INPUT` and `SEER_EXTENDED` to the two input files on your own machine. Optional `PNEN_PROJECT` must point to this repository root.

The original input is a tab-delimited, unformatted 30-column case-listing export with a header. Its column order is:

```
id, record, sequence, total, sex, race, year, marital, site, laterality,
hist, behavior, confirmation, source, grade_old, grade_new, stage,
surgery_old, surgery_new, no_surgery, radiation, chemo, months, days,
flag, vital, cancer, other, cod, age
```

The scripts decode year by adding 1800 and histology by adding 8000. Do not use this decoding on formatted or previously decoded values. The source export covers adult pulmonary neuroendocrine tumors diagnosed in 2004–2023; cohort selection is performed in code, including earliest-record selection before survival and study-window exclusions.

The supplementary input is a separate formatted-label 37-column export. Its first 30 fields correspond to the same records as the original export. The seven additional fields are harmonized tumor size, CS tumor size, AJCC sixth-edition stage group, T, N, M and schema. The first two fields identify patient and record for linkage. The script checks the original export size and the analyzed cohort size; adaptation to another export requires explicit review of those cohort checks. Patient identifiers are used locally for linkage only.

## Execution

```
Rscript check_code.R
Rscript run_analysis.R --list
Rscript run_analysis.R --core
Rscript run_analysis.R --all
```

`--core` runs the demographic OS, competing-risk, audit, clinical standardization and clinical mortality analyses. `--all` additionally runs supporting sensitivity, age/RMST, treatment/era and plotting analyses, in dependency order. Scripts run in separate R processes so temporary objects and random-number states do not leak between analyses. Full bootstrap runs can be computationally intensive. `check_code.R` parses all scripts and checks packages and internal script dependencies; it does not run an analysis or require data.

`work/pnen_adjusted.R` supplies definitions to other scripts and need not be run independently. `run_analysis.R` lists the exact execution order and stops at the first failure. It records package versions locally before execution.

## Analysis map

| Scripts | Purpose |
| --- | --- |
| `pnen_formal_OS.R`, `pnen_adjusted.R` | Demographic standardization, weighted Kaplan–Meier survival, joint bootstrap and stage-change contrasts |
| `pnen_formal_competing.R` | Weighted Aalen–Johansen mortality components |
| `pnen_formal_sensitivity.R`, `pnen_sensitivity_contrasts.R` | Cohort restrictions and simultaneous contrast families |
| `pnen_horizon1.R`, `pnen_horizon_event_check.R` | One-year prediction horizon and event support |
| `pnen_formal_descriptive.R`, `pnen_reporting_descriptives.R`, `pnen_formal_audit.R` | Cohort summaries and validation |
| `pnen_extension_audit.R`, `pnen_age_support_audit.R`, `pnen_extension_descriptive.R`, `pnen_extension_models.R` | Treatment composition, surgery strata and diagnosis-era extensions |
| `pnen_lung_additions.R` | Age-stratified probabilities and conditional restricted mean survival time |
| `pnen_value_audit.R` | Individual probability intervals, cause counts and fixed-reference support audit |
| `pnen_clinical_standardization.R` | C0/C1/C2 models adding AJCC stage and tumor size |
| `pnen_clinical_competing.R` | Clinical mortality components and common-bootstrap direct comparisons |
| `pnen_figures.R`, `pnen_extension_figures.R`, `pnen_lung_figures.R`, `pnen_tlcr_figure.R` | Statistical figures; the last script produces the final conditional-prognosis figure |

## Estimands and reproducibility

Landmarks are diagnosis and years 1, 2, 3 and 5; the principal horizon is the following 3 years. Differences are LCNEC minus SCLC. The `IPTW` label retained in code denotes inverse probability weighting for histology-group standardization, not a treatment-effect analysis. Overlap weighting targets a different population. Mortality states include survival, cancer-related death, other-cause death and unknown-cause death.

The demographic primary analyses use seed 20261005 and 2,000 bootstrap samples. Restricted-cohort analyses use seed 20261007 and 1,000 samples; treatment/era extensions use seed 20261010 and 1,000; the one-year horizon uses seed 20261013 and 1,000; age/RMST uses seed 20261014 and 1,000; clinical models and mortality comparisons use seed 20261017 and 1,000. Clinical mortality code checks that the OS bootstrap estimates reproduce those from the clinical survival models. Exploratory extensions retain separate uncertainty families and are not a prospectively registered confirmatory analysis.

## Generated files

Execution creates local files under `outputs/` and may create a patient-level intermediate `work/pnen_extension_data.rds`. None are part of this repository. The root-level allowlist in `.gitignore` excludes every non-code file, including generated estimates, figures and patient-level intermediates. Never force-add such files to a public repository.

Some internal output directories retain the historical name `pnen_manuscript`; they hold generated statistical outputs when code is run and do not contain a distributed manuscript. The original analysis scripts are preserved without changes to estimators, cohort rules, seeds or bootstrap counts.
