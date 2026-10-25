# Statistical research report — Social Trust & Wellbeing (GSS 2024)

**Opeyemi O. Agboola · 2024 General Social Survey Release 3a**

## Objective and data
Examine whether interpersonal trust is associated with three-category subjective happiness while analyzing institutional confidence separately. The verified original source has 3,986 cases: 3,309 core sample respondents and 677 AmeriSpeak supplement cases. The core sample is the analysis universe.

## Methods
Verified respondent eligibility; handled GSS reserved numeric codes; survey-weighted descriptive proportions; adjusted proportional-odds ordinal regression; 67 strata / 134 PSUs; WTSSNRPS weighting; stratified PSU score-based covariance; comparisons with alternative weights and interview-mode adjustment. Age, sex, degree and marital status were covariates. All estimates are observational, not causal.

## Core results
| Analysis quantity | Value |
| --- | --- |
| Valid trust responses | 946 |
| Joint trust/happiness respondents | 936 |
| Adjusted model complete cases | 911 |
| Weighted % reporting not too happy, cautious group | 27.4% |
| Weighted % reporting not too happy, high-trust group | 18.5% |
| Adjusted ordinal OR (cautious vs high trust) | 1.09 (95% CI 0.71–1.66) |
| Adjusted ordinal p | 0.698 |
| Joint exploratory trust test p | approximately 0.404 |

**Finding:** Weighted descriptive groups differed, but trust was not conclusively associated with happiness after adjustment.

## Reproducibility and limitations
Original GSS SPSS source: https://gss.norc.org/get-the-data/spss.html . Original data fingerprint SHA256: 048e0939f4c800fab15dd7bf71f9c7fb0ea0b4cac217c3b27d59ca280fdf34c6 . Numerical results and QA generated in Python. SPSS and R replication scripts prepared but native SPSS Viewer output and survey-model R replication have not been executed. GSS questionnaire routing, item eligibility, missingness, survey inference and proportional odds diagnostics constrain interpretation. No original respondent-level dataset has been published. Independent research study; not a client commission or an official NORC analysis.

Source: Davern et al. (2026), General Social Survey 1972–2024, Release 3a, NORC.