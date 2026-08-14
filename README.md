# MDCalc Entropy Analysis Code

This private repository contains only the executable analysis code and pinned Python dependencies for the MDCalc information-yield paper. Manuscripts, submission files, figures, source workbooks, cleaned data, and generated results are intentionally excluded.

**Code package:** `0.1.0`

**Full private working repository:** https://github.com/ShuhanCS/mdcalc-entropy

## Run the analysis

Use Python 3.12:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python analyze.py C:\path\to\mdcalc_entropy_source_v0.2.0.xlsx
```

The program verifies the expected source-workbook SHA-256 before analysis and writes generated data, tables, figures, and result summaries into the local repository. Those outputs are ignored by Git and must not be committed here.

## Repository boundary

Keep this repository code-only. Paper development, private data, manuscript drafts, citation files, and journal submission packages belong in the full private working repository linked above.
