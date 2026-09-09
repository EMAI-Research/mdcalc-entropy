# Zenodo data and code deposit

Publish the public reproducibility release for the submitted BMJ Evidence-Based Medicine manuscript, preserving the frozen analysis. User authorized the Zenodo deposit on September 8, 2026 (local time).

- Use the existing public repository as the release boundary. Include aggregate extracted data, audit data, code, result summaries, dictionary and scoped licenses. Exclude the source workbook, which contains third-party text, and submission correspondence.
- Confirm that the four public CSV datasets match the submitted medRxiv reproducibility archive. Keep CC BY 4.0 for first-party data/documentation and MIT for code.
- Update the title and deposit metadata to the 14-author submitted manuscript. Bump public package 0.2.0 to 0.3.0 for the archival release; scientific inputs and results are unchanged.
- Verify hashes, run the public reproduction commands in an isolated folder, build the ZIP, and publish it to Zenodo. Record only a DOI returned by a confirmed publication.
- Update the public repository and he-lab tracker with the DOI, or record the exact deposit blocker if Zenodo remains unavailable. Do not alter the already-submitted journal files or substitute the dataset DOI for the pending medRxiv DOI.
