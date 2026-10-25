# Social trust and subjective wellbeing: GSS 2024

Original quantitative analysis of the NORC 2024 General Social Survey, using the core sample (N = 3,309). The trust/happiness joint universe contained 936 respondents; 911 contributed complete information to adjusted models.

The survey-weighted descriptive contrast was 27.4% reporting "not too happy" among cautious respondents versus 18.5% among those trusting most people. Adjusted ordinal odds ratio: 1.09 (95% CI 0.71–1.66), inconclusive.

**Tools actually executed:** Python, NumPy, pandas, SciPy, statsmodels and survey-aware inference. SPSS and R replication scripts were prepared but not executed; no native SPSS Viewer export is claimed.

**Official source:** https://gss.norc.org/get-the-data/spss.html

**Limitations:** questionnaire routing; complete-case exclusions; proportional-odds assumptions; no causal inference.


## Published evidence

[Executed ordinal model](analysis/final_survey_models.py) | [SPSS reader](analysis/read_sav_fallback.py) | [Data audit](analysis/initial_audit.py) | [SPSS syntax](analysis/spss_descriptives.sps) | [Research report](RESEARCH_REPORT.md) | [Weighted estimates](results/weighted_descriptive_results.csv).

To reproduce, obtain the official 2024 SPSS GSS dataset and save it as data/raw/GSS2024.sav. The survey estimates were generated with Python; native SPSS Viewer and independent R-survey model results are not represented as completed.
