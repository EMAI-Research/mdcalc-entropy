# MDCalc entropy analysis

Public analysis code and de-identified aggregate data for the manuscript **“Meaning of clinical calculator results: a cross-sectional analysis of the MDCalc catalogue.”**

- Code release: `0.3.0`
- Journal submission: BMJ Evidence-Based Medicine, `bmjebm-2026-115053`
- Scientific analysis: frozen BMJ `0.14.3` analysis retained in EBM submission package `0.1.5`
- Data freeze: 11 August 2026
- Python: 3.12

The frozen analysis contains 482 evaluations of 407 calculators linked to a catalogue of 847 MDCalc calculators. Median reduction in starting uncertainty was 14.9%, or 0.093 bits.

## Run from the included derived data

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python build_result_support.py
python decision_relevance.py
python redesign_figures.py
python figure5_ruleout_support.py
```

These commands recalculate the support estimates and build the final figures from `data/mdcalc_evaluations_primary.csv`. Each script includes consistency checks and stops if the frozen counts do not match.

## Rebuild from the source workbook

`analyze.py` reconstructs the cleaned datasets, tables, figures, and result summaries from the frozen source workbook:

```powershell
python analyze.py C:\path\to\mdcalc_entropy_source_v0.2.0.xlsx
```

The program accepts only the workbook with SHA-256 `02348f2e95039be35943e135d27b12010ba382d589d03623a497baef35dbc514`. The workbook contains third-party source material and is not redistributed in this public repository. The included CSV files are the frozen, derived analytic datasets used in the manuscript.

## Contents

- `analyze.py`: deterministic reconstruction, analysis, checks, tables, and supporting figures.
- `redesign_figures.py`: final publication-figure generator.
- `build_result_support.py`: post-negative and post-positive risk estimates with Wilson 95% intervals.
- `decision_relevance.py` and `figure5_ruleout_support.py`: supporting decision-threshold analyses and figure code.
- `data/`: catalogue, primary, expanded, and row-level audit datasets.
- `results_summary.*`, `decision_relevance_results.json`, and `result_support_summary.json`: frozen result summaries.
- `DATA_DICTIONARY.md`: definitions for the derived data fields.
- `MANIFEST.json`: file sizes and SHA-256 hashes for the release.

The unit of analysis is one calculator-study-outcome-threshold evaluation. The files do not rank calculators, and the results are research outputs rather than clinical decision support. See `LICENSE.md` and `LICENSING_STATUS.md` for the scoped code and data licences.

## Archival deposit

Zenodo deposit metadata are in `.zenodo.json`; the deposit record is in `docs/zenodo-deposit.md`. No Zenodo DOI has been issued for this release yet. The dataset DOI will be separate from the medRxiv preprint DOI.

The four CSV datasets match the submitted reproducibility archive row for row. Code and data remain unchanged from public release 0.2.0; release 0.3.0 updates the manuscript title, author metadata and archival documentation.
