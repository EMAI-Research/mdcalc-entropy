# Zenodo deposit record

Status: prepared for upload; not deposited. No record ID or DOI has been issued.

The author authorized publication of the data and code archive. On September 9, 2026 UTC, Zenodo returned browser timeouts and HTTP 504; direct access also timed out. Publication remains pending access to the service and an authenticated account.

## Deposit contents

- Package: public reproducibility release 0.3.0.
- Title and 14 creators: `.zenodo.json`, copied from the submitted manuscript author list; author email addresses are excluded.
- Source: https://github.com/EMAI-Research/mdcalc-entropy.
- Primary data: 482 evaluations; expanded data: 489; audit data: 494; catalogue: 847.
- All four CSVs match the submitted medRxiv reproducibility archive row for row.
- Source workbook excluded because it includes third-party text. Original workbook SHA-256: `02348f2e95039be35943e135d27b12010ba382d589d03623a497baef35dbc514`.
- First-party data/documentation: CC BY 4.0. Original code: MIT. These are aggregate research data, not individual patient data.

## Validation

All four documented public reproduction commands passed in an isolated Python 3.12 environment with the pinned requirements: `build_result_support.py`, `decision_relevance.py`, `redesign_figures.py`, and `figure5_ruleout_support.py`. Creator metadata includes all 14 submitted authors and the 13 supplied ORCIDs; no author email addresses are included. The archive manifest verifies every packaged file except the manifest itself. Scientific code and data match public release 0.2.0.

## Complete the deposit

Upload `mdcalc-entropy-data-code-v0.3.0.zip` at https://zenodo.org/uploads/new using `.zenodo.json` metadata. Verify the file checksum and all 14 creators before publishing. Record the published record URL and DOI here, add the DOI to the README and he-lab tracker, and provide the data citation to BMJ through its available update mechanism. Do not use the dataset DOI in the journal's preprint DOI field.
