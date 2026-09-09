# Raw data (not included)

This folder is intentionally empty of survey microdata. `scripts/07_empirical_main.R`
expects one file here:

```
raw/ZA7600_v3-0-0.dta
```

## How to obtain it

The empirical analysis (H1, H3, and the auxiliary channel test) uses the
**ISSP 2019 Social Inequality V** module (GESIS study **ZA7600**), version
3.0.0 or later, Stata (`.dta`) format.

1. Register for a free GESIS account and request the dataset at
   https://doi.org/10.4232/1.14009 (or search "ZA7600" at
   https://www.gesis.org/en/issp/data-and-documentation/social-inequality/2019).
   Citation: ISSP Research Group (2022). International Social Survey
   Programme: Social Inequality V -- ISSP 2019. GESIS, Cologne. ZA7600
   Data file Version 3.0.0, https://doi.org/10.4232/1.14009.
2. Download the Stata-format file and place it at `raw/ZA7600_v3-0-0.dta`.
3. Run `scripts/07_empirical_main.R`.

## Why the data itself is not included here

GESIS's usage regulations prohibit redistributing licensed microdata to
third parties (see
https://www.gesis.org/fileadmin/user_upload/Usage_regulations.pdf). This
repository therefore ships only the code that operates on the data, plus
derived, respondent-anonymous outputs (model coefficients, country-level
aggregates, simulation output) in `data/` and `results/` — nothing that
reconstructs individual ISSP responses.

All simulation-based results (Figures 1–4, S1–S4 and the corresponding
data files) are entirely self-contained and require no external data.
